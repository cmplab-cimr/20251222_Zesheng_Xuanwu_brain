# -*- coding: utf-8 -*-
"""器官 × 化学类 富集热图 (Figure 4 Panel C / 设计方案 v2)

行 = 化学类, 列 = 器官, colorbar 在右侧竖直。
格色 = 富集比 (蓝耗竭-白1.0-红富集), cell 文字 = FDR 星号。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...figure_plotting_package.basic_shape_elements.shapes import Rectangle
from ...figure_plotting_package.common.core_plotting_functions import HeatmapValueFormat
from ...common.common_figure_materials import rdbu_diverging_cmap, make_colorbar

HeatmapDataFigure = GeneralElements.HeatmapDataFigure


class OrganClassEnrichmentHeatmap(CompositeFigure):
    """器官 × 化学类 富集热图 (转置: 行=类, 列=器官, cbar 在右)

    figure_data_parameter_dict:
        organs: list[str]
        classes: list[str]        — 前 lipid_n 个脂质, 其余非脂质
        ratio_matrix: 2D array    — [organ, class] 富集比 (内部转置)
        fdr_matrix: 2D array      — [organ, class] FDR
        lipid_n: int              — 脂质数 (水平分隔线)
    """

    margin_left = 0.22
    margin_bottom = 0.08
    margin_right = 0.08
    margin_top = 0.04
    cbar_w = 0.03
    gap_cbar = 0.03

    def __init__(self, figure_data_parameter_dict, **kwargs):
        organs = default_parameter_extract(
            figure_data_parameter_dict, 'organs', [], pop=True)
        classes = default_parameter_extract(
            figure_data_parameter_dict, 'classes', [], pop=True)
        ratio_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'ratio_matrix', np.array([[1.0]]), pop=True)
        fdr_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'fdr_matrix', np.array([[1.0]]), pop=True)
        lipid_n = default_parameter_extract(
            figure_data_parameter_dict, 'lipid_n', 0, pop=True)
        if not organs or not classes:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        # 转置: 行=类, 列=器官
        ratio_matrix = np.array(ratio_matrix, dtype=float).T
        fdr_matrix = np.array(fdr_matrix, dtype=float).T
        n_cls, n_org = ratio_matrix.shape

        dev = float(np.nanmax(np.abs(ratio_matrix - 1.0))) if ratio_matrix.size else 0.5
        dev = max(dev, 0.3)
        data_lim = (1.0 - dev, 1.0 + dev)

        from ...common.common_figure_materials import fdr_star
        fdr_arr = fdr_matrix
        def _star_format(x, y, decimal_value):
            return fdr_star(fdr_arr[y, x]) or None
        data_value_text_format_dict = {
            ParameterName.basic_number_format_str: _star_format,
            'font_size': 7,
        }

        # 脂质/非脂质水平分隔线
        highlight_config = [(Rectangle, {
            'center': Vector((n_org - 1) / 2.0, lipid_n - 0.5),
            'width': n_org,
            'height': 0.03,
            'face_color': '#555555',
            'edge_width': 0,
            'z_order': 5,
        })]

        heatmap_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: rdbu_diverging_cmap(),
            ParameterName.data_nested_list: ratio_matrix,
            ParameterName.data_lim_pair: data_lim,
            ParameterName.figure_config_dict: {
                ParameterName.data_value_text_format_dict: data_value_text_format_dict,
                ParameterName.x_tick_label_format_dict: {
                    ParameterName.font_size: 8,
                },
                ParameterName.y_tick_label_format_dict: {
                    ParameterName.font_size: 8,
                },
            },
            ParameterName.hidden_axis_boundary: False,
            ParameterName.x_tick_labels_list: [o.capitalize()[:7] for o in organs],
            ParameterName.y_tick_labels_list: classes,
            ParameterName.highlight: True,
            ParameterName.highlight_config: highlight_config,
        }
        plot_w = max(0.4, n_org * 0.10 + 0.10)
        plot_h = max(0.4, n_cls * 0.085 + 0.08)
        heatmap_obj = HeatmapDataFigure(
            heatmap_param,
            bottom_left=Vector(self.margin_left, self.margin_bottom),
            size=Vector(plot_w, plot_h))

        # 右侧竖直 colorbar
        cbar_bl = Vector(self.margin_left + plot_w + self.gap_cbar, self.margin_bottom)
        cbar_sz = Vector(self.cbar_w, plot_h)
        cmap = rdbu_diverging_cmap()
        cbar_obj = make_colorbar(cmap, data_lim[0], data_lim[1], cbar_bl, cbar_sz,
                                label='Enrichment', orientation='vertical')

        subfigure_element_dict = {
            'heatmap': {'heatmap': heatmap_obj},
            'colorbar': {'enrich_cbar': cbar_obj},
        }
        total_width = self.margin_left + plot_w + self.gap_cbar + self.cbar_w + self.margin_right
        total_height = self.margin_bottom + plot_h + self.margin_top
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
