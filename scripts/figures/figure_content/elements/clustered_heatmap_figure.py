# -*- coding: utf-8 -*-
"""差异代谢物聚类热图 Panel B

布局（顶→底，左→右，参考论文样式）：
  顶部: 列 dendrogram (scipy) + 组别色条
  中心: 主热图 (z-score, 灰→白→青发散)
  左侧: 代谢物分类色条 (Class I, 强调几类脂质)
  右侧: 强调脂质连线 + 标签
  底部: 横向 colorbar

架构: ClusteredHeatmapFigure(CompositeFigure) 只组装 subfigure_element_dict，
      不重写 draw（复用框架递归渲染）。各区块做成 DataFigureAxes 子类
      （复用 generate_mpl_axes 拿到正确 figure rect），在绝对 figure-width
      坐标空间定位，与 CompositeFigure 的 total_width/total_height 同空间，
      由 move_and_scale 统一缩放（参考 oplsda_grid_figure.py 模式）。
"""

from scipy.cluster import hierarchy
import matplotlib.colors as mcolors

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, np,
)
from ...figure_plotting_package.basic_shape_elements.composite_figure_and_axes import DataFigureAxes
from ...common.common_figure_materials import gray_teal_diverging_cmap, make_colorbar, CommonColorDict as CCD

HeatmapDataFigure = GeneralElements.HeatmapDataFigure


# 脂质类配色/判断 (lipid_class_color, is_lipid_class) 统一定义在 _colors.py, Panel B/C 共用


# ============================================================
# DataFigureAxes 子类: 每个区块一个独立 axes
# ============================================================

class DendrogramAxes(DataFigureAxes):
    """顶部列聚类 dendrogram, 水平方向覆盖样本列范围。"""

    def __init__(self, col_linkage, bottom_left, size, **kwargs):
        self.col_linkage = col_linkage
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        hierarchy.dendrogram(
            self.col_linkage, orientation='top', no_labels=True,
            ax=ax, color_threshold=0,
            link_color_func=lambda *a: '#333333',
            above_threshold_color='#333333')
        # 不设 xlim: scipy dendrogram 叶在 x=5+10*i (自动 xlim [0,10n]), 均匀分布,
        # 相对位置 (i+0.5)/n 与热图 cell (xlim -0.5..n-0.5) 一致 → figure 上天然对齐
        for coll in ax.collections:    # dendrogram 线条变细
            coll.set_linewidth(0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax, transform


class AnnotationBar(DataFigureAxes):
    """通用注释色条。

    orientation='top':  一行 n_cols 个方块（组别色条, 按列）
    orientation='left': 一列 n_rows 个方块（Class I 色条, 按行）
    """

    def __init__(self, labels, color_map, orientation, bottom_left, size, **kwargs):
        self.labels = labels
        self.color_map = color_map
        self.orientation = orientation
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        n = len(self.labels)
        if n == 0:
            return ax, transform
        colors = [self.color_map.get(label, CCD.light_gray) for label in self.labels]
        cmap_obj = mcolors.ListedColormap(colors)
        if self.orientation == 'top':
            data = np.arange(n).reshape(1, n)
            ax.imshow(data, aspect='auto', cmap=cmap_obj,
                      vmin=-0.5, vmax=n - 0.5, interpolation='nearest')
            ax.set_xlim(-0.5, n - 0.5)        # 按列对齐
            ax.set_ylim(-0.5, 0.5)
        else:  # 'left'
            data = np.arange(n).reshape(n, 1)
            ax.imshow(data, aspect='auto', cmap=cmap_obj,
                      vmin=-0.5, vmax=n - 0.5, interpolation='nearest')
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(n - 0.5, -0.5)        # 按行对齐 (与热图 y 翻转一致, 行 0 在顶)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax, transform


class ColorBarAxes(DataFigureAxes):
    """底部横向 z-score colorbar (自画 imshow, 带数值标签)。"""

    def __init__(self, cmap, vmin, vmax, bottom_left, size, label='Z-score', **kwargs):
        self.cmap = cmap
        self.vmin = vmin
        self.vmax = vmax
        self.label = label
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        n = 256
        gradient = np.linspace(self.vmin, self.vmax, n).reshape(1, n)
        ax.imshow(gradient, aspect='auto', cmap=self.cmap,
                  vmin=self.vmin, vmax=self.vmax, interpolation='nearest')
        tick_pos = np.linspace(0, n - 1, 5)
        tick_val = np.linspace(self.vmin, self.vmax, 5)
        ax.set_xticks(tick_pos)
        ax.set_xticklabels([f'{v:.1f}' for v in tick_val], fontsize=5)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color('#cccccc')
        return ax, transform


class HighlightLabelAxes(DataFigureAxes):
    """右侧强调脂质: 从热图行中点拉折线到右侧标签。"""

    def __init__(self, highlight_indices, highlight_labels, n_rows, bottom_left, size, **kwargs):
        self.indices = highlight_indices
        self.labels = highlight_labels
        self.n_rows = n_rows
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        ax.set_xlim(0, 1)                          # 逻辑坐标: 0=热图右边缘
        ax.set_ylim(self.n_rows - 0.5, -0.5)       # 必须镜像热图 y 翻转
        for row_idx, label in zip(self.indices, self.labels):
            y = row_idx                           # cell center (热图行坐标)
            ax.plot([0, 0.55], [y, y], color='#333333', lw=0.6, zorder=5)
            ax.text(0.58, y, label, fontsize=5, va='center', ha='left',
                    fontfamily='DejaVu Sans', zorder=6, clip_on=False)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax, transform


# ============================================================
# ClusteredHeatmapFigure: 组装各区块
# ============================================================

class ClusteredHeatmapFigure(CompositeFigure):
    """差异代谢物聚类热图 (dendrogram + 注释条 + 主热图 + 强调标签 + colorbar)

    figure_data_parameter_dict:
        z_matrix_ordered:  np.ndarray       — (n_rows, n_cols) 已按行列聚类排序的 z-score
        col_linkage:       np.ndarray       — 列聚类 scipy linkage (画顶部 dendrogram)
        groups_ordered:    list[str]        — 每列分组 (按列聚类顺序)
        group_colors:      dict[str,str]    — {group: hex}
        metabolite_labels: list[str]        — 行标签 (主热图 y tick)
        sample_labels:     list[str]        — 列标签 (主热图 x tick, 可空)
        class_i_ordered:   list[str]        — 每行 Class I 标签 (按行聚类顺序)
        class_color_map:   dict[str,str]    — {class_label: hex}
        highlight_indices: list[int]        — 强调代谢物在 z_matrix_ordered 中的行下标
        highlight_labels:  list[str]        — 强调代谢物显示名
        z_lim:             tuple            — (vmin, vmax), 默认 (-2, 2)
        cmap:              colormap         — 默认 gray_teal_diverging_cmap
    """

    # 布局常量 (figure-width fraction, 与 total_width 同空间)
    left_class_bar_w = 0.018
    gap_class_heat = 0.01      # 左侧分类色条与主图的横向 gap
    right_label_w = 0.18
    top_dendro_h = 0.11        # 顶部 dendrogram 高度
    top_group_h = 0.012
    bottom_cbar_h = 0.03
    bottom_margin = 0.01       # colorbar 距 figure 底
    gap_dendro_group = 0.004
    gap_group_heat = 0.004
    gap_heat_cbar = 0.03       # 主图与 colorbar 的纵向 gap

    def __init__(self, figure_data_parameter_dict, **kwargs):
        z_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'z_matrix_ordered', np.array([[0]]), pop=True)
        col_linkage = default_parameter_extract(
            figure_data_parameter_dict, 'col_linkage', None, pop=True)
        groups_ordered = default_parameter_extract(
            figure_data_parameter_dict, 'groups_ordered', [], pop=True)
        group_colors = default_parameter_extract(
            figure_data_parameter_dict, 'group_colors', {}, pop=True)
        metabolite_labels = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_labels', [], pop=True)
        sample_labels = default_parameter_extract(
            figure_data_parameter_dict, 'sample_labels', [], pop=True)
        class_i_ordered = default_parameter_extract(
            figure_data_parameter_dict, 'class_i_ordered', [], pop=True)
        class_color_map = default_parameter_extract(
            figure_data_parameter_dict, 'class_color_map', {}, pop=True)
        highlight_indices = default_parameter_extract(
            figure_data_parameter_dict, 'highlight_indices', [], pop=True)
        highlight_labels = default_parameter_extract(
            figure_data_parameter_dict, 'highlight_labels', [], pop=True)
        z_lim = default_parameter_extract(
            figure_data_parameter_dict, 'z_lim', (-2, 2), pop=True)
        cmap = default_parameter_extract(
            figure_data_parameter_dict, 'cmap', gray_teal_diverging_cmap(), pop=True)

        n_rows, n_cols = z_matrix.shape
        if n_rows == 0 or n_cols == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        # 各区块绝对尺寸 (figure-width fraction)
        heatmap_w = max(0.45, min(0.75, 0.02 * n_cols + 0.3))
        heatmap_h = max(0.3, min(0.95, 0.012 * n_rows + 0.15))

        left_pad = self.left_class_bar_w + self.gap_class_heat
        right_pad = self.right_label_w
        top_pad = (self.top_dendro_h + self.gap_dendro_group
                   + self.top_group_h + self.gap_group_heat)
        bottom_pad = self.bottom_margin + self.bottom_cbar_h + self.gap_heat_cbar

        total_width = left_pad + heatmap_w + right_pad
        total_height = top_pad + heatmap_h + bottom_pad

        # 中心热图 (原点左下)
        hm_bl = Vector(left_pad, bottom_pad)
        hm_sz = Vector(heatmap_w, heatmap_h)
        # 顶部按列区块: x 与热图一致 → 像素按列对齐
        col_x = left_pad
        col_w = heatmap_w
        grp_y = bottom_pad + heatmap_h + self.gap_group_heat
        grp_bl = Vector(col_x, grp_y)
        grp_sz = Vector(col_w, self.top_group_h)
        dnd_y = grp_y + self.top_group_h + self.gap_dendro_group
        dnd_bl = Vector(col_x, dnd_y)
        dnd_sz = Vector(col_w, self.top_dendro_h)
        # 左侧分类色条: x=[0, left_class_bar_w], 与主图留 gap_class_heat; y 与热图一致
        cls_bl = Vector(0, bottom_pad)
        cls_sz = Vector(self.left_class_bar_w, heatmap_h)
        # 底部 colorbar: 距 figure 底 bottom_margin, 与主图留 gap_heat_cbar
        cbar_bl = Vector(col_x, self.bottom_margin)
        cbar_sz = Vector(col_w, self.bottom_cbar_h)
        # 右侧强调区: y 与热图一致
        rlbl_bl = Vector(left_pad + heatmap_w, bottom_pad)
        rlbl_sz = Vector(self.right_label_w, heatmap_h)

        # 主热图 (框架 HeatmapDataFigure, 自动 xlim/ylim 对齐 + tick label)
        heatmap_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: cmap,
            ParameterName.data_nested_list: z_matrix,
            ParameterName.data_lim_pair: z_lim,
            ParameterName.figure_config_dict: {},
            ParameterName.x_tick_labels_list: None,   # 不显示下方小鼠编号
            ParameterName.y_tick_labels_list: None,   # 不显示左侧代谢物名
            ParameterName.hidden_axis_boundary: False,
            ParameterName.hidden_axis_tick: True,     # 隐藏 tick 标记
        }
        heatmap_obj = HeatmapDataFigure(heatmap_param, bottom_left=hm_bl, size=hm_sz)

        # 组装 subfigure_element_dict (dict 保序 = 渲染 z-order, highlight 须在 heatmap 之后)
        subfigure_element_dict = {}
        if col_linkage is not None:
            subfigure_element_dict['dendrogram'] = {
                'col_dendro': DendrogramAxes(col_linkage, dnd_bl, dnd_sz)}
        subfigure_element_dict['annotation'] = {
            'group_bar': AnnotationBar(groups_ordered, group_colors, 'top', grp_bl, grp_sz),
            'class_bar': AnnotationBar(class_i_ordered, class_color_map, 'left', cls_bl, cls_sz),
        }
        subfigure_element_dict['heatmap'] = {'main_heat': heatmap_obj}
        if highlight_indices:
            subfigure_element_dict['highlight'] = {
                'right_lbl': HighlightLabelAxes(
                    highlight_indices, highlight_labels, n_rows, rlbl_bl, rlbl_sz)}
        subfigure_element_dict['colorbar'] = {
            'z_cbar': make_colorbar(cmap, z_lim[0], z_lim[1], cbar_bl, cbar_sz, label='Z-score')}

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
