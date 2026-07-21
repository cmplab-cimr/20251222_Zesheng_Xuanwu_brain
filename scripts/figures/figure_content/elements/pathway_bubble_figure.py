# -*- coding: utf-8 -*-
"""跨器官通路富集点图

Panel D: organs × pathways 点图
点大小 = pathway impact, 颜色 = 方向 / 显著性
通路按「显著组织数」从大到小排序, 通路名旋转 90 度避免重叠

用 DataFigureAxes 子类自画 scatter（参考 Panel B 模式），
完全控制 tick label rotation 与排序，不依赖框架 tick label 旋转。
"""

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, np,
)
from ...figure_plotting_package.basic_shape_elements.composite_figure_and_axes import DataFigureAxes
from ...common.common_figure_materials import CommonColorDict as CCD


class PathwayBubbleAxes(DataFigureAxes):
    """通路气泡散点 axes（matplotlib 原生 scatter + 旋转 90 的 x tick label）。"""

    def __init__(self, x, y, sizes, colors, n_pathway, n_tissue, pathways, tissues,
                 bottom_left, size, **kwargs):
        self.x, self.y, self.sizes, self.colors = x, y, sizes, colors
        self.n_pathway, self.n_tissue = n_pathway, n_tissue
        self.pathways, self.tissues = pathways, tissues
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        if len(self.x):
            ax.scatter(self.x, self.y, s=self.sizes, c=self.colors, alpha=0.8,
                       edgecolors='#555555', linewidths=0.4, zorder=3)
        ax.set_xticks(range(self.n_pathway))
        ax.set_xticklabels(
            [p[:24] + '..' if len(p) > 25 else p for p in self.pathways],
            rotation=90, fontsize=4.5, ha='center', fontfamily='DejaVu Sans')
        ax.set_yticks(range(self.n_tissue))
        ax.set_yticklabels(
            [t.capitalize() for t in self.tissues],
            fontsize=6, fontfamily='DejaVu Sans')
        ax.set_xlim(-0.5, self.n_pathway - 0.5)
        ax.set_ylim(-0.5, self.n_tissue - 0.5)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
        return ax, transform


class PathwayBubbleFigure(CompositeFigure):
    """跨器官通路富集点图

    figure_data_parameter_dict:
        pathway_df: DataFrame  — columns: tissue, pathway_name, p_value, pathway_impact, median_log2FC, ...
        max_pathways: int  — 最多显示的通路数
        color_by: str  — 'direction' 或 'significance'
        tissue_order: list[str]  — 组织排列
    """

    # 旋转 90 的通路标签需要较大底部 margin
    margin_left = 0.12
    margin_bottom = 0.28
    margin_right = 0.04
    margin_top = 0.04

    def __init__(self, figure_data_parameter_dict, **kwargs):
        pathway_df = default_parameter_extract(
            figure_data_parameter_dict, 'pathway_df', None, pop=True)
        max_pathways = default_parameter_extract(
            figure_data_parameter_dict, 'max_pathways', 20, pop=True)
        color_by = default_parameter_extract(
            figure_data_parameter_dict, 'color_by', 'direction', pop=True)
        tissue_order = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_order', None, pop=True)

        if pathway_df is None or pathway_df.empty:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        # 通路排序：按「显著组织数 (p<0.05)」从大到小, 平局按 min p 升序
        stats = (pathway_df.assign(_sig=pathway_df['p_value'] < 0.05)
                 .groupby('pathway_name')
                 .agg(_sig_count=('_sig', 'sum'), _min_p=('p_value', 'min'))
                 .sort_values(['_sig_count', '_min_p'], ascending=[False, True]))
        pathways = stats.head(max_pathways).index.tolist()

        tissues = tissue_order or sorted(pathway_df['tissue'].unique())
        tissues = [t for t in tissues if t in pathway_df['tissue'].values]
        df = pathway_df[pathway_df['pathway_name'].isin(pathways)]
        pathways = [p for p in pathways if p in set(df['pathway_name'])]  # 保持显著度排序

        n_tissue = len(tissues)
        n_pathway = len(pathways)
        if n_tissue == 0 or n_pathway == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        tissue_y = {t: i for i, t in enumerate(tissues)}
        pathway_x = {p: i for i, p in enumerate(pathways)}

        xs, ys, sizes, colors = [], [], [], []
        for _, row in df.iterrows():
            if row['tissue'] in tissue_y and row['pathway_name'] in pathway_x:
                xs.append(pathway_x[row['pathway_name']])
                ys.append(tissue_y[row['tissue']])
                impact = row.get('pathway_impact', row.get('impact', 0.1))
                sizes.append(max(15, min(180, float(impact) * 700)))
                if color_by == 'direction':
                    colors.append(CCD.teal if row.get('median_log2FC', 0) > 0 else CCD.gray)
                else:
                    colors.append(CCD.teal if -np.log10(max(row['p_value'], 1e-300)) > 1.3 else CCD.light_gray)

        plot_width = max(0.5, n_pathway * 0.055)
        plot_height = max(0.3, n_tissue * 0.06)
        total_width = self.margin_left + plot_width + self.margin_right
        total_height = self.margin_bottom + plot_height + self.margin_top

        axes_obj = PathwayBubbleAxes(
            np.array(xs) if xs else np.array([]),
            np.array(ys) if ys else np.array([]),
            np.array(sizes) if sizes else np.array([]),
            colors, n_pathway, n_tissue, pathways, tissues,
            bottom_left=Vector(self.margin_left, self.margin_bottom),
            size=Vector(plot_width, plot_height))

        subfigure_element_dict = {'bubble': {'bubble': axes_obj}}
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
