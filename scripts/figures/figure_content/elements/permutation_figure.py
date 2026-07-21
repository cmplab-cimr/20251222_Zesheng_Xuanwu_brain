# -*- coding: utf-8 -*-
"""OPLS-DA 置换检验网格 (Figure S4B 模型验证)

用框架 ScatterDataFigure 画 (风格统一), 10 组织小网格, 每格:
  permuted R2Y(x) vs permuted Q2(y) 灰散点 + 真模型青点 + 回归线 + p 值标题。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, LineStyle,
)
from ...figure_plotting_package.common.color import ColorConfig
from ...common.common_figure_materials import CommonColorDict as CCD

ScatterDataFigure = GeneralElements.ScatterDataFigure

TISSUE_DISPLAY = {
    'cortex': 'Cortex', 'hippocampus': 'Hippocampus', 'hypothalamus': 'Hypothalamus',
    'liver': 'Liver', 'kidney': 'Kidney', 'fat': 'Adipose',
    'gastrocnemius': 'Gastrocnemius', 'quadriceps': 'Quadriceps', 'ileum': 'Ileum',
    'serum': 'Serum',
}


def _build_perm_scatter(perm_r2y, perm_q2, orig_r2y, orig_q2, p_r2y, p_q2, tissue,
                        bottom_left, size):
    n_perm = len(perm_r2y)
    x_all = np.concatenate([perm_r2y, [orig_r2y]])
    y_all = np.concatenate([perm_q2, [orig_q2]])
    marker_size = np.concatenate([np.full(n_perm, 6.0), [20.0]])
    marker_color = [CCD.gray] * n_perm + [CCD.teal]
    x_lim_pad = (np.max(x_all) - np.min(x_all)) * 0.1
    y_lim_pad = (np.max(y_all) - np.min(y_all)) * 0.1
    x_lim = Vector(float(np.min(x_all) - x_lim_pad), float(np.max(x_all) + x_lim_pad))
    y_lim = Vector(float(np.min(y_all) - y_lim_pad), float(np.max(y_all) + y_lim_pad))

    data_dict = {
        ParameterName.x_value_array: x_all,
        ParameterName.y_value_array: (y_all, None),
        ParameterName.marker_size: marker_size,
        ParameterName.marker_color: marker_color,
        ParameterName.scatter_param_dict: {
            ParameterName.alpha: 0.6,
            ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order,
        },
    }
    # 回归线
    line_config = {
        **DataFigureConfig.common_line_param_dict_generator(),
        ParameterName.edge_color: ColorConfig.gray,
        ParameterName.z_order: DataFigureConfig.line_z_order,
        ParameterName.edge_style: LineStyle.normal,
    }
    scatter_line = []
    if n_perm >= 2 and np.ptp(perm_r2y) > 0:
        z = np.polyfit(perm_r2y, perm_q2, 1)
        xs = np.linspace(x_lim[0], x_lim[1], 20)
        scatter_line.append([xs, np.polyval(z, xs), dict(line_config)])

    title = f'{TISSUE_DISPLAY.get(tissue, tissue)}\npR2Y={p_r2y:.3f} pQ2={p_q2:.3f}'
    scatter_param = {
        ParameterName.ax_bottom_left_list: [Vector(0, 0)],
        ParameterName.ax_size_list: [Vector(1, 1)],
        ParameterName.color_dict: None,
        ParameterName.data_nested_list: [data_dict],
        ParameterName.scatter_line: [scatter_line],
        ParameterName.figure_config_dict: {
            ParameterName.subplot_name_text_format_dict: {
                **DataFigureConfig.common_subplot_text_format_dict_generator(),
                ParameterName.font_size: 6,
            },
        },
        ParameterName.subplot_name_list: [title],
        ParameterName.text_axis_loc_pair: Vector(0.5, 1.12),
        ParameterName.x_lim_list: [x_lim],
        ParameterName.x_label_list: ['perm R2Y'],
        ParameterName.x_ticks_list: [[]],
        ParameterName.y_lim_list: [y_lim],
        ParameterName.y_label_list: ['perm Q2'],
        ParameterName.y_ticks_list: [[]],
        ParameterName.legend: False,
    }
    return ScatterDataFigure(scatter_param, bottom_left=bottom_left, size=size)


class PermutationGridFigure(CompositeFigure):
    """多组织置换检验网格

    figure_data_parameter_dict:
        perm_results: dict  — {tissue: {perm_r2y, perm_q2, original_r2y, original_q2, p_r2y, p_q2}}
        tissue_order: list[str]
    """

    subplot_width = 0.22
    subplot_height = 0.22
    subplot_gap_x = 0.03
    subplot_gap_y = 0.05

    def __init__(self, figure_data_parameter_dict, **kwargs):
        perm_results = default_parameter_extract(
            figure_data_parameter_dict, 'perm_results', {}, pop=True)
        tissue_order = default_parameter_extract(
            figure_data_parameter_dict, 'tissue_order', list(perm_results.keys()), pop=True)
        valid = [t for t in tissue_order if t in perm_results]
        n = len(valid)
        if n == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        ncols = min(5, n)
        nrows = (n + ncols - 1) // ncols
        total_width = ncols * self.subplot_width + (ncols - 1) * self.subplot_gap_x
        total_height = nrows * self.subplot_height + (nrows - 1) * self.subplot_gap_y

        subfigure_element_dict = {}
        for idx, t in enumerate(valid):
            pr = perm_results[t]
            row, col = divmod(idx, ncols)
            x0 = col * (self.subplot_width + self.subplot_gap_x)
            y0 = (nrows - 1 - row) * (self.subplot_height + self.subplot_gap_y)
            scatter_obj = _build_perm_scatter(
                np.asarray(pr['perm_r2y']), np.asarray(pr['perm_q2']),
                pr['original_r2y'], pr['original_q2'],
                pr.get('p_r2y', 0), pr.get('p_q2', 0), t,
                bottom_left=Vector(x0, y0),
                size=Vector(self.subplot_width, self.subplot_height))
            name = f'perm_{t}'
            subfigure_element_dict[name] = {name: scatter_obj}

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
