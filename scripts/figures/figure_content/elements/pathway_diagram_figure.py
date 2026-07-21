# -*- coding: utf-8 -*-
"""磷脂重塑 / 溶血磷脂—类花生酸 通路示意图 (Figure 4 Panel F1 / 设计方案"改法一")

PC → PLA2 切割 → LPC/LPE/LPS/LPA + 游离 PUFA(AA/DHA) → COX/LOX → 氧化脂质
节点按脑区(海马) log2FC 方向上色 (红=上调, 蓝=下调, 灰=酶/无数据)。
描述性着色, 非富集检验, 无偏倚问题。用 DataFigureAxes 自画。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...figure_plotting_package.basic_shape_elements.composite_figure_and_axes import DataFigureAxes

# 方向色 (红=上调, 蓝=下调, 灰=酶/无数据/不显著)
UP_COLOR = '#d85a30'
DOWN_COLOR = '#3a7ca5'
NEUTRAL_COLOR = '#CCCCCC'
INSIG_COLOR = '#EEEEEE'


def _node_color(log2fc, is_enzyme=False, fc_threshold=0.1):
    if is_enzyme or log2fc is None or not np.isfinite(log2fc):
        return NEUTRAL_COLOR
    if abs(log2fc) < fc_threshold:
        return INSIG_COLOR
    return UP_COLOR if log2fc > 0 else DOWN_COLOR


class PathwayDiagramAxes(DataFigureAxes):
    """通路示意图: 节点(色=海马方向) + 箭头 + 文字标签。"""

    def __init__(self, nodes, edges, bottom_left, size, **kwargs):
        self.nodes = nodes      # [{label, x, y, log2fc, is_enzyme}]
        self.edges = edges      # [(from_idx, to_idx)]
        super().__init__(bottom_left, size, **kwargs)

    def draw(self, fig, parent_ax=None, parent_transformation=None):
        ax, transform = self.generate_mpl_axes(fig, parent_ax)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        # 箭头 (先画, 在节点下)
        for fi, ti in self.edges:
            fn, tn = self.nodes[fi], self.nodes[ti]
            ax.annotate(
                '', xy=(tn['x'], tn['y']), xytext=(fn['x'], fn['y']),
                arrowprops=dict(arrowstyle='->', color='#999999', lw=0.8,
                                shrinkA=8, shrinkB=8),
                zorder=1)
        # 节点
        for node in self.nodes:
            col = _node_color(node.get('log2fc'), node.get('is_enzyme', False))
            ax.scatter([node['x']], [node['y']], s=350, c=col,
                       edgecolors='#555555', linewidths=0.5, zorder=3)
            ax.text(node['x'], node['y'], node['label'],
                    fontsize=3.5, ha='center', va='center',
                    fontfamily='DejaVu Sans', zorder=4, wrap=True)
        # 图例 (方向色)
        for i, (lab, col) in enumerate([('Up', UP_COLOR), ('Down', DOWN_COLOR), ('n.s.', INSIG_COLOR)]):
            ax.scatter([0.02], [0.06 - i * 0.04], s=60, c=col, edgecolors='#555', linewidths=0.3, zorder=5)
            ax.text(0.04, 0.06 - i * 0.04, lab, fontsize=4, va='center', ha='left', fontfamily='DejaVu Sans', zorder=5)
        ax.set_axis_off()
        return ax, transform


class PathwayDiagramFigure(CompositeFigure):
    """磷脂重塑通路示意图

    figure_data_parameter_dict:
        nodes: list[dict]  — {label, x, y, log2fc, is_enzyme}
        edges: list[tuple] — [(from_idx, to_idx)]
    """

    margin = 0.04

    def __init__(self, figure_data_parameter_dict, **kwargs):
        nodes = default_parameter_extract(
            figure_data_parameter_dict, 'nodes', [], pop=True)
        edges = default_parameter_extract(
            figure_data_parameter_dict, 'edges', [], pop=True)
        if not nodes:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        total_width = 1.0
        total_height = 0.85
        axes_obj = PathwayDiagramAxes(
            nodes, edges,
            bottom_left=Vector(self.margin, self.margin),
            size=Vector(total_width - 2 * self.margin, total_height - 2 * self.margin))
        subfigure_element_dict = {'pathway': {'pathway': axes_obj}}
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
