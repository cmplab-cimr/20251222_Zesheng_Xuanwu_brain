# -*- coding: utf-8 -*-
"""关键脂质代谢物跨器官柱图 + vs 表型相关

Panel F: 上部跨组织分组柱状图 + 下部代谢物×表型相关热图
- 嵌套 CrossTissueBarFigure (上部) + PhenotypeCorrelationFigure (下部)
- CompositeFigure 内部用 move_and_scale 定位两个子 CompositeFigure
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure,
    default_parameter_extract,
)
from .cross_tissue_bar_figure import CrossTissueBarFigure
from .phenotype_corr_figure import PhenotypeCorrelationFigure


class LipidCrossTissuePhenotypeFigure(CompositeFigure):
    """关键脂质代谢物跨器官柱图 + vs 表型相关组合图。

    figure_data_parameter_dict:
        # —— CrossTissueBarFigure 参数 ——
        metabolites: list[str]          — 代谢物 ID 列表
        values: dict                    — {met_id: {tissue: {group: array}}}
        metabolite_names: dict          — {met_id: display_name}
        tissue_order: list[str]         — 组织排列顺序
        group_colors: dict              — {group: color_hex}
        group_names: tuple              — (group1, group2)

        # —— PhenotypeCorrelationFigure 参数 ——
        r_matrix: np.ndarray           — (n_mets, n_phenotypes) 相关系数矩阵
        metabolite_labels_corr: list[str] — 相关热图行标签
        phenotype_labels: list[str]    — 表型名列表
    """

    GAP = 0.05           # 上下两部分间距 (figure-width fraction)
    TOP_FRAC = 0.55      # 上部柱图占总高度的比例

    def __init__(self, figure_data_parameter_dict, **kwargs):
        # —— 提取参数 ——
        metabolites = default_parameter_extract(
            figure_data_parameter_dict, 'metabolites', [], pop=True)
        values = default_parameter_extract(
            figure_data_parameter_dict, 'values', {}, pop=True)
        metabolite_names = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_names', {}, pop=True)
        tissue_order = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_order', [], pop=True)
        group_colors = default_parameter_extract(
            figure_data_parameter_dict, 'group_colors', {}, pop=True)
        group_names = default_parameter_extract(
            figure_data_parameter_dict, 'group_names', ('AD', 'Sema'), pop=True)

        r_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'r_matrix', np.array([[0.]]), pop=True)
        metabolite_labels_corr = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_labels_corr', [], pop=True)
        phenotype_labels = default_parameter_extract(
            figure_data_parameter_dict, 'phenotype_labels', [], pop=True)

        n_mets = len(metabolites)
        has_bar = n_mets > 0 and len(values) > 0
        has_corr = r_matrix.size > 0 and r_matrix.shape[0] > 0

        if not has_bar and not has_corr:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        # —— 创建两个子 CompositeFigure ——
        bar_fig = None
        corr_fig = None
        bar_w = bar_h = 0.0
        corr_w = corr_h = 0.0

        if has_bar:
            bar_params = {
                'metabolites': metabolites,
                'values': values,
                'metabolite_names': metabolite_names,
                'tissue_order': tissue_order,
                'group_colors': group_colors,
                'group_names': group_names,
            }
            bar_fig = CrossTissueBarFigure(figure_data_parameter_dict=bar_params)
            bar_w, bar_h = bar_fig.size.x, bar_fig.size.y

        if has_corr:
            r_matrix = np.asarray(r_matrix, dtype=float)
            corr_params = {
                'r_matrix': r_matrix,
                'metabolite_labels': metabolite_labels_corr,
                'phenotype_labels': phenotype_labels,
                'tissue': 'hippocampus',
            }
            corr_fig = PhenotypeCorrelationFigure(figure_data_parameter_dict=corr_params)
            corr_w, corr_h = corr_fig.size.x, corr_fig.size.y

        # —— 计算总尺寸 ——
        total_width = max(bar_w, corr_w)
        total_height = bar_h + (self.GAP if has_bar and has_corr else 0) + corr_h

        # —— 定位子元素 ——
        # corr_fig 在底部 (原点), bar_fig 在 corr 上方
        if corr_fig is not None:
            # corr_fig 中心对齐
            if corr_w < total_width:
                corr_fig.move_and_scale(
                    scale=1, bottom_left_offset=Vector((total_width - corr_w) / 2, 0))

        if bar_fig is not None:
            bar_y = corr_h + (self.GAP if has_corr else 0)
            if bar_w < total_width:
                bar_fig.move_and_scale(
                    scale=1, bottom_left_offset=Vector((total_width - bar_w) / 2, bar_y))
            else:
                bar_fig.move_and_scale(scale=1, bottom_left_offset=Vector(0, bar_y))

        # —— 组装 element_dict ——
        element_dict = {}
        if bar_fig is not None:
            element_dict['bar_section'] = {'cross_tissue_bar': bar_fig}
        if corr_fig is not None:
            element_dict['corr_section'] = {'phenotype_corr': corr_fig}

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
