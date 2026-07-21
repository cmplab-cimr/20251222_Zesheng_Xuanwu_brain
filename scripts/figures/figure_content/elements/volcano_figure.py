# -*- coding: utf-8 -*-
"""火山图 Panel C

继承共享代码 VolcanoScatterDataFigure 保持轴/阈值线风格统一, 仅替换 marker_size/marker_color:
- 显著脂质 = CCD.teal (与 Panel B 脂质色一致)
- 显著非脂质 = CCD.gray
- 不显著 = CCD.light_gray
- 所有点同样大小; 仅 HEATMAP_HIGHLIGHT_METABOLITES 里的重要脂质连线标名
"""

import numpy as np

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements, merge_axis_format_dict, ColorConfig, LineStyle,
)
from ...figure_elements.metabolism_volcano_figure import VolcanoScatterDataFigure
from ...common.common_figure_materials import CommonColorDict as CCD, is_lipid_class

ScatterDataFigure = GeneralElements.ScatterDataFigure


class LipidVolcanoScatterDataFigure(VolcanoScatterDataFigure):
    """火山图: 显著脂质 CCD.teal / 显著非脂质 CCD.gray / 不显著 CCD.light_gray;
    统一点大小; 仅重要脂质标名。继承 VolcanoScatterDataFigure 风格统一。"""

    def __init__(self, figure_data_parameter_dict, bottom_left, size, **kwargs):
        ax_total_bottom_left = Vector(0, 0)
        ax_total_size = Vector(1, 1) - ax_total_bottom_left

        (log_fold, log_p, log_fold_thr, log_p_thr) = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_data, '', pop=True, force=True)
        lipid_mask = default_parameter_extract(
            figure_data_parameter_dict, 'lipid_mask', None, pop=True)
        highlight_indices = default_parameter_extract(
            figure_data_parameter_dict, 'highlight_indices', [], pop=True)
        names = default_parameter_extract(
            figure_data_parameter_dict, 'names', [], pop=True)
        marker_size = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.marker_size, 10, pop=True)
        x_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.common_x_label, 'log2(Fold Change)', pop=True)
        y_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.common_y_label, '-log10(p)', pop=True)
        new_figure_config_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_config_dict, {}, pop=True)

        # 显著 = p 超阈值(不分上下调); marker_color per-point
        significant = log_p > log_p_thr
        if lipid_mask is None:
            lipid_mask = [False] * len(log_fold)
        marker_color_list = [
            CCD.teal if (significant[i] and lipid_mask[i])
            else (CCD.gray if significant[i] else CCD.light_gray)
            for i in range(len(log_fold))
        ]

        x_min = np.min(log_fold)
        x_max = np.max(log_fold)
        y_max = np.max(log_p)
        x_lim = Vector(x_min, x_max) + Vector(-1, 1) * 0.1 * (x_max - x_min)
        y_lim = (-0.9, 1.1 * y_max)

        # supplementary_text: 仅重要脂质
        supplementary_text_list = []
        supplementary_text_loc_list = []
        x_span = x_lim[1] - x_lim[0]
        y_span = y_lim[1] - y_lim[0]
        for i in highlight_indices:
            if 0 <= i < len(names):
                supplementary_text_list.append(names[i])
                supplementary_text_loc_list.append(np.array([
                    (log_fold[i] - x_lim[0]) / x_span,
                    (log_p[i] - y_lim[0]) / y_span + 0.02,
                ]))

        current_data_dict = {
            ParameterName.x_value_array: log_fold,
            ParameterName.y_value_array: (log_p, None),
            ParameterName.marker_size: marker_size,
            ParameterName.marker_color: marker_color_list,
            ParameterName.scatter_param_dict: {
                ParameterName.alpha: 0.8,
                ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order,
            },
        }

        supplementary_text_format_dict = {
            **DataFigureConfig.common_supplementary_text_config_dict,
            ParameterName.font_size: 4.5,
            ParameterName.width: 0.1,
            ParameterName.height: 0.05,
        }
        flux_title_text_format_dict = {
            **DataFigureConfig.common_subplot_text_format_dict_generator(),
            ParameterName.font_size: 7,
        }
        figure_config_dict = {
            ParameterName.y_label_format_dict: merge_axis_format_dict(
                {}, {ParameterName.axis_label_distance: 0.04},
                new_figure_config_dict, ParameterName.y_label_format_dict),
            ParameterName.x_tick_label_format_dict: merge_axis_format_dict(
                {}, {ParameterName.axis_tick_label_distance: 0.008},
                new_figure_config_dict, ParameterName.x_tick_label_format_dict),
            ParameterName.x_label_format_dict: merge_axis_format_dict(
                {}, {ParameterName.axis_label_distance: 0.025},
                new_figure_config_dict, ParameterName.x_label_format_dict),
            ParameterName.subplot_name_text_format_dict: flux_title_text_format_dict,
            ParameterName.line_param_dict: {
                **DataFigureConfig.common_line_param_dict_generator(),
                ParameterName.z_order: DataFigureConfig.line_z_order,
                ParameterName.edge_style: LineStyle.dash,
            },
            ParameterName.supplementary_text_format_dict: supplementary_text_format_dict,
        }
        common_dash_line_config_dict = {
            **DataFigureConfig.common_line_param_dict_generator(),
            ParameterName.edge_color: ColorConfig.gray,
            ParameterName.z_order: DataFigureConfig.line_z_order,
            ParameterName.edge_style: LineStyle.dash,
        }
        scatter_line_pair_list = [[
            [x_lim[0], x_lim[1]], [log_p_thr, log_p_thr], dict(common_dash_line_config_dict),
        ], [
            [log_fold_thr, log_fold_thr], [y_lim[0], y_lim[1]], dict(common_dash_line_config_dict),
        ], [
            [-log_fold_thr, -log_fold_thr], [y_lim[0], y_lim[1]], dict(common_dash_line_config_dict),
        ]]

        figure_data_parameter_dict = {
            ParameterName.ax_bottom_left_list: [ax_total_bottom_left],
            ParameterName.ax_size_list: [ax_total_size],
            ParameterName.color_dict: None,
            ParameterName.data_nested_list: [current_data_dict],
            ParameterName.figure_config_dict: figure_config_dict,
            ParameterName.scatter_line: [scatter_line_pair_list],
            ParameterName.x_lim_list: [x_lim],
            ParameterName.x_ticks_list: [[]],
            ParameterName.x_label_list: [x_label],
            ParameterName.y_lim_list: [y_lim],
            ParameterName.y_ticks_list: [[]],
            ParameterName.y_label_list: [y_label],
            ParameterName.supplementary_text_list: [supplementary_text_list],
            ParameterName.supplementary_text_loc_list: [supplementary_text_loc_list],
            ParameterName.legend: False,
            **figure_data_parameter_dict,
        }
        # 跳过 VolcanoScatterDataFigure.__init__ (其 marker_color 用 color_dict), 直接调 ScatterDataFigure
        ScatterDataFigure.__init__(self, figure_data_parameter_dict, bottom_left, size, **kwargs)


class MetabolomicsVolcanoFigure(CompositeFigure):
    """火山图 — 显著脂质 CCD.teal/显著非脂质 CCD.gray/不显著 CCD.light_gray, 统一点大小, 重要脂质标名

    figure_data_parameter_dict:
        log2fc, neg_log10p: np.ndarray
        classes: list — 每个代谢物的 Class I 标签
        names: list[str]
        highlight_indices: list[int] — 重要脂质在数组中的下标
        thresholds: tuple — (fc_threshold, p_threshold)
        marker_size: float — 统一点大小
    """

    def __init__(self, figure_data_parameter_dict, **kwargs):
        log2fc = default_parameter_extract(
            figure_data_parameter_dict, 'log2fc', np.array([]), pop=True)
        neg_log10p = default_parameter_extract(
            figure_data_parameter_dict, 'neg_log10p', np.array([]), pop=True)
        classes = default_parameter_extract(
            figure_data_parameter_dict, 'classes', [], pop=True)
        names = default_parameter_extract(
            figure_data_parameter_dict, 'names', [], pop=True)
        highlight_indices = default_parameter_extract(
            figure_data_parameter_dict, 'highlight_indices', [], pop=True)
        fc_threshold, p_threshold = default_parameter_extract(
            figure_data_parameter_dict, 'thresholds', (1.5, 0.05), pop=True)
        marker_size = default_parameter_extract(
            figure_data_parameter_dict, 'marker_size', 10, pop=True)

        lipid_mask = [is_lipid_class(c) for c in classes]
        log_fc_thr = np.log2(fc_threshold) if fc_threshold > 1 else fc_threshold
        log_p_thr = -np.log10(p_threshold) if 0 < p_threshold < 1 else p_threshold

        total_width = 1.0
        total_height = 0.9
        scatter_obj = LipidVolcanoScatterDataFigure(
            {
                ParameterName.figure_data: (log2fc, neg_log10p, log_fc_thr, log_p_thr),
                'lipid_mask': lipid_mask,
                'highlight_indices': highlight_indices,
                'names': names,
                ParameterName.marker_size: marker_size,
                ParameterName.common_x_label: 'log2(Fold Change)',
                ParameterName.common_y_label: '-log10(p)',
            },
            bottom_left=Vector(0, 0),
            size=Vector(total_width, total_height),
        )
        subfigure_element_dict = {'volcano': {'volcano': scatter_obj}}
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
