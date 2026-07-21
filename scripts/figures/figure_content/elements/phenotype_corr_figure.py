# -*- coding: utf-8 -*-
"""代谢物–表型相关热图

Panel G: 代谢物 (y) × 表型 (x), 颜色 = r, 标记显著性
"""

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
    np,
)

HeatmapDataFigure = GeneralElements.HeatmapDataFigure

from ...common.common_figure_materials import rdbu_diverging_cmap, make_colorbar


class PhenotypeCorrelationFigure(CompositeFigure):
    """代谢物–表型相关热图

    figure_data_parameter_dict:
        r_matrix: np.ndarray  — 相关系数矩阵 (n_metabolites × n_phenotypes)
        p_matrix: np.ndarray  — p-value 矩阵
        metabolite_labels: list[str]  — 行标签
        phenotype_labels: list[str]  — 列标签
        tissue: str  — 组织名
    """

    title_height = 0.04
    cbar_h = 0.03
    gap_cbar = 0.02

    def __init__(self, figure_data_parameter_dict, **kwargs):
        r_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'r_matrix', np.array([[0]]), pop=True)
        metabolite_labels = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_labels', [], pop=True)
        phenotype_labels = default_parameter_extract(
            figure_data_parameter_dict, 'phenotype_labels', [], pop=True)
        tissue = default_parameter_extract(
            figure_data_parameter_dict, 'tissue', '', pop=True)

        n_rows, n_cols = r_matrix.shape
        if n_rows == 0 or n_cols == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False)
            return

        # 热图尺寸
        heatmap_width = max(0.3, n_cols * 0.12)
        heatmap_height = max(0.3, min(0.8, n_rows * 0.03))

        # 标签缩写
        met_labels = [
            l[:20] + '..' if len(str(l)) > 20 else str(l)
            for l in metabolite_labels
        ]

        heatmap_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: rdbu_diverging_cmap(),
            ParameterName.data_nested_list: r_matrix,
            ParameterName.data_lim_pair: (-1, 1),
            ParameterName.figure_config_dict: {},
            ParameterName.x_tick_labels_list: phenotype_labels,
            ParameterName.y_tick_labels_list: met_labels,
            ParameterName.hidden_axis_boundary: False,
        }

        heatmap_obj = HeatmapDataFigure(
            heatmap_param,
            bottom_left=Vector(0, self.cbar_h + self.gap_cbar),
            size=Vector(heatmap_width, heatmap_height),
        )

        # 底部 colorbar
        cbar_obj = make_colorbar(rdbu_diverging_cmap(), -1, 1, Vector(0, 0),
                                Vector(heatmap_width, self.cbar_h), label='r')

        subfigure_element_dict = {
            'heatmap': {'heatmap': heatmap_obj},
            'colorbar': {'r_cbar': cbar_obj},
        }

        total_width = heatmap_width + 0.02
        total_height = self.cbar_h + self.gap_cbar + heatmap_height

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
