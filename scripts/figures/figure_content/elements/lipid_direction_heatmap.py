# -*- coding: utf-8 -*-
"""已注释脂质代谢物跨器官方向热图

Panel D (D1): metabolite × organ 的 log2FC 方向热图
- 左侧: 代谢物名称 (行标签)
- 中心: 方向热图 (蓝=下调, 白=不变, 红=上调)
- 底部: 横向 colorbar (log2FC)

布局架构: LipidMetaboliteDirectionFigure(CompositeFigure) 组装 subfigure_element_dict，
不重写 draw。
"""

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance

from ...common.config import (
    ParameterName, Vector, CompositeFigure,
    default_parameter_extract, GeneralElements,
)
from ...common.common_figure_materials import rdbu_diverging_cmap, make_colorbar

HeatmapDataFigure = GeneralElements.HeatmapDataFigure


class LipidMetaboliteDirectionFigure(CompositeFigure):
    """已注释脂质代谢物跨器官方向热图。

    figure_data_parameter_dict:
        direction_matrix: np.ndarray     — (n_metabolites, n_tissues) log2FC 值
        metabolite_labels: list[str]     — 行标签 (代谢物显示名, 左侧显示)
        tissue_labels: list[str]         — 列标签 (组织缩写)
        cmap: colormap | None           — 默认 rdbu_diverging_cmap
        data_lim: tuple | None          — (vmin, vmax)，默认 auto
        cluster_rows: bool              — 默认 True，行聚类排序
    """

    # 布局常量 (figure-width fraction)
    LABEL_WIDTH = 0.22               # 左侧代谢物名预留宽度
    GAP_LABEL_HEAT = 0.01
    BOTTOM_CBAR_H = 0.030
    BOTTOM_MARGIN = 0.005
    GAP_HEAT_CBAR = 0.025
    TOP_MARGIN = 0.02

    def __init__(self, figure_data_parameter_dict, **kwargs):
        direction_matrix = default_parameter_extract(
            figure_data_parameter_dict, 'direction_matrix', np.array([[0.]]), pop=True)
        metabolite_labels = default_parameter_extract(
            figure_data_parameter_dict, 'metabolite_labels', [], pop=True)
        tissue_labels = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_labels', [], pop=True)
        cmap = default_parameter_extract(
            figure_data_parameter_dict, 'cmap', rdbu_diverging_cmap(), pop=True)
        data_lim = default_parameter_extract(
            figure_data_parameter_dict, 'data_lim', None, pop=True)
        cluster_rows = default_parameter_extract(
            figure_data_parameter_dict, 'cluster_rows', True, pop=True)

        direction_matrix = np.asarray(direction_matrix, dtype=float)
        n_rows, n_cols = direction_matrix.shape

        if n_rows == 0 or n_cols == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        # —— 行聚类 ——
        if cluster_rows and n_rows > 1:
            dmat = np.nan_to_num(direction_matrix, nan=0.0)
            row_linkage = hierarchy.linkage(distance.pdist(dmat, metric='euclidean'), method='ward')
            row_order = [int(x) for x in hierarchy.leaves_list(row_linkage)]
            direction_matrix = direction_matrix[row_order, :]
            metabolite_labels = [metabolite_labels[i] for i in row_order]

        # —— 色标范围 ——
        if data_lim is None:
            finite = direction_matrix[np.isfinite(direction_matrix)]
            abs_finite = np.abs(finite)
            vmax = float(np.percentile(abs_finite, 95)) if abs_finite.size else 1.0
            vmax = max(vmax, 0.5)
            vmax = min(vmax, 10.0)
            data_lim = (-vmax, vmax)

        # —— 各区块尺寸 ——
        heatmap_w = max(0.35, min(0.65, 0.06 * n_cols + 0.25))
        heatmap_h = max(0.25, min(0.90, 0.025 * n_rows + 0.15))

        left_pad = self.LABEL_WIDTH + self.GAP_LABEL_HEAT
        bottom_pad = self.BOTTOM_MARGIN + self.BOTTOM_CBAR_H + self.GAP_HEAT_CBAR

        total_width = left_pad + heatmap_w
        total_height = self.TOP_MARGIN + heatmap_h + bottom_pad

        # 中心热图
        hm_bl = Vector(left_pad, bottom_pad)
        hm_sz = Vector(heatmap_w, heatmap_h)
        # 底部 colorbar
        cbar_bl = Vector(left_pad, self.BOTTOM_MARGIN)
        cbar_sz = Vector(heatmap_w, self.BOTTOM_CBAR_H)

        # —— 主热图 ——
        heatmap_param = {
            ParameterName.ax_bottom_left_list: Vector(0, 0),
            ParameterName.ax_size_list: Vector(1, 1),
            ParameterName.heatmap_cmap: cmap,
            ParameterName.data_nested_list: direction_matrix,
            ParameterName.data_lim_pair: data_lim,
            ParameterName.figure_config_dict: {},
            ParameterName.x_tick_labels_list: tissue_labels,
            ParameterName.y_tick_labels_list: metabolite_labels,  # 左侧代谢物名
            ParameterName.hidden_axis_boundary: False,
        }
        heatmap_obj = HeatmapDataFigure(heatmap_param, bottom_left=hm_bl, size=hm_sz)

        # —— 组装 ——
        subfigure_element_dict = {
            'heatmap': {'direction_heat': heatmap_obj},
            'colorbar': {
                'log2fc_cbar': make_colorbar(cmap, data_lim[0], data_lim[1], cbar_bl, cbar_sz, label='log2FC'),
            },
        }

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
