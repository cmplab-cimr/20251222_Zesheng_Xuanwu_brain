# -*- coding: utf-8 -*-
"""代谢物 × 表型 相关散点网格 (Figure 4 Panel G / Figure S4)

每小图: x=代谢物丰度, y=表型值, 按 AD(gray)/Sema(teal) 分组上色 + 全体线性回归线。
用框架 ScatterDataFigure (像 CrossTissueBarFigure 一次传多个小图数据)。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, LineStyle,
)
from ...common.common_figure_materials import CommonColorDict as CCD

ScatterDataFigure = GeneralElements.ScatterDataFigure


class CorrelationScatterGridFigure(CompositeFigure):
    """相关散点网格: 多个小图 (ScatterDataFigure), 每图散点+回归线+标题。

    figure_data_parameter_dict:
        panels: list[dict] — 每图 {x, y, group_colors, fit_x, fit_y, title, x_label, y_label}
            x: np.array — 代谢物丰度 (全体动物)
            y: np.array — 表型值
            group_colors: list[str] — per-point 颜色 (CCD.gray / CCD.teal)
            fit_x/fit_y: np.array — 回归线 x/y
            title: str — 子图标题 (代谢物名 + r/p)
    """

    subplot_width = 0.30
    subplot_height = 0.30
    subplot_gap_x = 0.04
    subplot_gap_y = 0.06

    def __init__(self, figure_data_parameter_dict, **kwargs):
        panels = default_parameter_extract(
            figure_data_parameter_dict, 'panels', [], pop=True)
        n = len(panels)
        if n == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        ncols = min(2, n)
        nrows = (n + ncols - 1) // ncols
        total_width = ncols * self.subplot_width + (ncols - 1) * self.subplot_gap_x
        total_height = nrows * self.subplot_height + (nrows - 1) * self.subplot_gap_y

        line_config = {
            **DataFigureConfig.common_line_param_dict_generator(),
            ParameterName.edge_color: '#333333',
            ParameterName.z_order: DataFigureConfig.line_z_order,
            ParameterName.edge_style: LineStyle.normal,
        }

        subfigure_element_dict = {}
        for idx, panel in enumerate(panels):
            row, col = divmod(idx, ncols)
            x0 = col * (self.subplot_width + self.subplot_gap_x)
            y0 = (nrows - 1 - row) * (self.subplot_height + self.subplot_gap_y)

            data_dict = {
                ParameterName.x_value_array: np.asarray(panel['x'], dtype=float),
                ParameterName.y_value_array: (np.asarray(panel['y'], dtype=float), None),
                ParameterName.marker_size: 12,
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
                        ParameterName.font_size: 6,
                    },
                },
                ParameterName.subplot_name_list: [panel.get('title', '')],
                ParameterName.text_axis_loc_pair: Vector(0.5, 1.08),
                ParameterName.x_label_list: [panel.get('x_label', '')],
                ParameterName.y_label_list: [panel.get('y_label', '')],
                ParameterName.x_ticks_list: [[]],
                ParameterName.y_ticks_list: [[]],
                ParameterName.legend: False,
            }
            obj = ScatterDataFigure(
                scatter_param,
                bottom_left=Vector(x0, y0),
                size=Vector(self.subplot_width, self.subplot_height))
            name = f'corr_{idx}'
            subfigure_element_dict[name] = {name: obj}

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
