# -*- coding: utf-8 -*-
"""认知–脂质相关热图 + 侧边 CCS 散点 (Figure 4 Panel G)

左: 代谢物(行) × 认知表型(列) Spearman r 热图 (rdbu 发散色), cell 星号 = raw p;
右: 头条代谢物 × CCS 散点 (AD=gray/Sema=teal + 回归线) + 竖直 colorbar。
热图模式照搬 organ_class_enrichment_heatmap, 散点照搬 correlation_scatter_grid_figure。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, LineStyle,
)
from ...common.common_figure_materials import rdbu_diverging_cmap, make_colorbar, fdr_star

HeatmapDataFigure = GeneralElements.HeatmapDataFigure
ScatterDataFigure = GeneralElements.ScatterDataFigure


class CognitionCorrelationFigure(CompositeFigure):
    """认知–脂质相关热图 + 侧边 CCS 散点。

    figure_data_parameter_dict:
        r_matrix: np.ndarray        — (n_metabolites × n_phenotypes) Spearman r
        p_matrix: np.ndarray        — 同形状 raw p_value (用于星号标注)
        metabolite_labels: list[str]
        phenotype_labels: list[str]
        inset_panel: dict | None    — 头条散点 (CorrelationScatterGridFigure 单 panel 格式:
                                      {x, y, group_colors, fit_x, fit_y, title, x_label, y_label})
        tissue: str
    """

    margin_left = 0.26      # 代谢物名空间
    margin_bottom = 0.10    # 表型 x 标签空间
    margin_right = 0.04
    margin_top = 0.04
    cell_w = 0.14
    cell_h = 0.055
    cbar_w = 0.025
    gap_cbar = 0.03
    gap_scatter = 0.06
    scatter_w = 0.32

    def __init__(self, figure_data_parameter_dict, **kwargs):
        r_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'r_matrix', np.array([[0.0]]), pop=True)
        p_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'p_matrix', np.ones_like(r_matrix), pop=True)
        metabolite_labels = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_labels', [], pop=True)
        phenotype_labels = default_parameter_extract(
            figure_data_parameter_dict, 'phenotype_labels', [], pop=True)
        inset_panel = default_parameter_extract(
            figure_data_parameter_dict, 'inset_panel', None, pop=True)
        _tissue = default_parameter_extract(figure_data_parameter_dict, 'tissue', '', pop=True)

        r_matrix = np.array(r_matrix, dtype=float)
        p_matrix = np.array(p_matrix, dtype=float)
        n_rows, n_cols = r_matrix.shape
        if n_rows == 0 or n_cols == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        met_labels = [str(l)[:22] + '..' if len(str(l)) > 22 else str(l) for l in metabolite_labels]

        # 星号 cell 文本 (raw p; n=12 pilot 下 BH-FDR 后多不显著, 图注声明 nominal p)
        def _star_format(x, y, decimal_value):
            return fdr_star(p_matrix[y, x]) or None

        data_value_text_format_dict = {
            ParameterName.basic_number_format_str: _star_format,
            'font_size': 8,
        }

        plot_w = max(0.35, n_cols * self.cell_w)
        plot_h = max(0.35, min(0.75, n_rows * self.cell_h))
        heatmap_bl = Vector(self.margin_left, self.margin_bottom)

        heatmap_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: rdbu_diverging_cmap(),
            ParameterName.data_nested_list: r_matrix,
            ParameterName.data_lim_pair: (-1, 1),
            ParameterName.figure_config_dict: {
                ParameterName.data_value_text_format_dict: data_value_text_format_dict,
                ParameterName.x_tick_label_format_dict: {ParameterName.font_size: 8},
                ParameterName.y_tick_label_format_dict: {ParameterName.font_size: 7},
            },
            ParameterName.hidden_axis_boundary: False,
            ParameterName.x_tick_labels_list: phenotype_labels,
            ParameterName.y_tick_labels_list: met_labels,
        }
        heatmap_obj = HeatmapDataFigure(
            heatmap_param, bottom_left=heatmap_bl, size=Vector(plot_w, plot_h))

        # 竖直 colorbar (热图右侧)
        cbar_bl = Vector(self.margin_left + plot_w + 0.01, self.margin_bottom)
        cbar_obj = make_colorbar(rdbu_diverging_cmap(), -1, 1, cbar_bl,
                                 Vector(self.cbar_w, plot_h), label='r', orientation='vertical')

        subfigure_element_dict = {
            'heatmap': {'heatmap': heatmap_obj},
            'colorbar': {'r_cbar': cbar_obj},
        }

        right_x = self.margin_left + plot_w + 0.01 + self.cbar_w + self.gap_scatter
        if inset_panel is not None:
            scatter_h = min(self.scatter_w, plot_h)
            scatter_bl = Vector(right_x, self.margin_bottom + (plot_h - scatter_h) / 2)
            subfigure_element_dict['scatter'] = {
                'cog_scatter': self._build_scatter(inset_panel, scatter_bl,
                                                   Vector(self.scatter_w, scatter_h))}
            total_width = right_x + self.scatter_w + self.margin_right
        else:
            total_width = right_x + self.margin_right
        total_height = self.margin_bottom + plot_h + self.margin_top

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)

    @staticmethod
    def _build_scatter(panel, bottom_left, size):
        """单散点 axes (复用 CorrelationScatterGridFigure 的 panel 格式)。"""
        line_config = {
            **DataFigureConfig.common_line_param_dict_generator(),
            ParameterName.edge_color: '#333333',
            ParameterName.z_order: DataFigureConfig.line_z_order,
            ParameterName.edge_style: LineStyle.normal,
        }
        data_dict = {
            ParameterName.x_value_array: np.asarray(panel['x'], dtype=float),
            ParameterName.y_value_array: (np.asarray(panel['y'], dtype=float), None),
            ParameterName.marker_size: 16,
            ParameterName.marker_color: panel['group_colors'],
            ParameterName.scatter_param_dict: {
                ParameterName.alpha: 0.8,
                ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order,
            },
        }
        scatter_line = []
        if 'fit_x' in panel and len(panel['fit_x']) >= 2:
            scatter_line.append([panel['fit_x'], panel['fit_y'], dict(line_config)])
        scatter_param = {
            ParameterName.ax_bottom_left_list: [Vector(0, 0)],
            ParameterName.ax_size_list: [Vector(1, 1)],
            ParameterName.color_dict: None,
            ParameterName.data_nested_list: [data_dict],
            ParameterName.scatter_line: [scatter_line],
            ParameterName.figure_config_dict: {
                ParameterName.subplot_name_text_format_dict: {
                    **DataFigureConfig.common_subplot_text_format_dict_generator(),
                    ParameterName.font_size: 7,
                },
            },
            ParameterName.subplot_name_list: [panel.get('title', '')],
            ParameterName.text_axis_loc_pair: Vector(0.5, 1.10),
            ParameterName.x_label_list: [panel.get('x_label', '')],
            ParameterName.y_label_list: [panel.get('y_label', '')],
            ParameterName.x_ticks_list: [[]],
            ParameterName.y_ticks_list: [[]],
            ParameterName.legend: False,
        }
        return ScatterDataFigure(scatter_param, bottom_left=bottom_left, size=size)
