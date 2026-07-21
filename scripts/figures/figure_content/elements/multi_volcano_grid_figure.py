# -*- coding: utf-8 -*-
"""多器官火山图网格 (Figure S4D 其余器官广度)

用共享 metabolism_volcano_figure.VolcanoScatterDataFigure (风格统一),
下丘脑/回肠/肾/股四头肌等组织火山图小网格, 每格阈值线 + 组织名标题。
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
)
from ...figure_elements.metabolism_volcano_figure import VolcanoScatterDataFigure
from ...common.common_figure_materials import CommonColorDict as CCD

TISSUE_DISPLAY = {
    'cortex': 'Cortex', 'hippocampus': 'Hippocampus', 'hypothalamus': 'Hypothalamus',
    'liver': 'Liver', 'kidney': 'Kidney', 'fat': 'Adipose',
    'gastrocnemius': 'Gastrocnemius', 'quadriceps': 'Quadriceps', 'ileum': 'Ileum',
    'serum': 'Serum',
}


class MultiVolcanoGridFigure(CompositeFigure):
    """多器官火山图小网格 (复用 VolcanoScatterDataFigure)

    figure_data_parameter_dict:
        volcano_data: dict  — {tissue: {log2fc, neg_log10p}}
        tissues: list[str]
        thresholds: tuple   — (fc_threshold, p_threshold)
    """

    subplot_width = 0.26
    subplot_height = 0.26
    subplot_gap_x = 0.03
    subplot_gap_y = 0.05

    def __init__(self, figure_data_parameter_dict, **kwargs):
        volcano_data = default_parameter_extract(
            figure_data_parameter_dict, 'volcano_data', {}, pop=True)
        tissues = default_parameter_extract(
            figure_data_parameter_dict, 'tissues', [], pop=True)
        fc_threshold, p_threshold = default_parameter_extract(
            figure_data_parameter_dict, 'thresholds', (1.5, 0.05), pop=True)

        valid = [t for t in tissues if t in volcano_data]
        n = len(valid)
        if n == 0:
            super().__init__({}, Vector(0, 0), Vector(0.1, 0.1), background=False, **kwargs)
            return

        fc_thr = np.log2(fc_threshold) if fc_threshold > 1 else fc_threshold
        p_thr = -np.log10(p_threshold) if 0 < p_threshold < 1 else p_threshold
        color_dict = {
            ParameterName.positive: CCD.teal,
            ParameterName.negative: CCD.gray,
            ParameterName.medium_data: CCD.light_gray,
        }

        ncols = min(4, n)
        nrows = (n + ncols - 1) // ncols
        total_width = ncols * self.subplot_width + (ncols - 1) * self.subplot_gap_x
        total_height = nrows * self.subplot_height + (nrows - 1) * self.subplot_gap_y

        subfigure_element_dict = {}
        for idx, t in enumerate(valid):
            d = volcano_data[t]
            row, col = divmod(idx, ncols)
            x0 = col * (self.subplot_width + self.subplot_gap_x)
            y0 = (nrows - 1 - row) * (self.subplot_height + self.subplot_gap_y)
            inner_param = {
                ParameterName.figure_data: (np.asarray(d['log2fc']), np.asarray(d['neg_log10p']),
                                            fc_thr, p_thr),
                ParameterName.color_dict: color_dict,
                ParameterName.common_x_label: 'log2(FC)',
                ParameterName.common_y_label: '-log10(p)',
                ParameterName.marker_size: 6,
                ParameterName.figure_config_dict: {
                    ParameterName.subplot_name_text_format_dict: {
                        **DataFigureConfig.common_subplot_text_format_dict_generator(),
                        ParameterName.font_size: 6,
                    },
                },
                ParameterName.subplot_name_list: [TISSUE_DISPLAY.get(t, t)],
                ParameterName.text_axis_loc_pair: Vector(0.5, 1.08),
            }
            vol_obj = VolcanoScatterDataFigure(
                inner_param,
                bottom_left=Vector(x0, y0),
                size=Vector(self.subplot_width, self.subplot_height))
            name = f'volcano_{t}'
            subfigure_element_dict[name] = {name: vol_obj}

        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
