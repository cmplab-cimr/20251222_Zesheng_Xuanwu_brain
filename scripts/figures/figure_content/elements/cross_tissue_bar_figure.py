# -*- coding: utf-8 -*-
"""关键代谢物跨组织柱状图

CrossTissueBarDataFigure(BarDataFigure): 一次性传入所有代谢物数据
CrossTissueBarFigure(CompositeFigure): 容器包裹
"""

from ...common.config import (
    ParameterName, Vector, CompositeFigure, DataFigureConfig,
    default_parameter_extract, GeneralElements,
    np,
)

BarDataFigure = GeneralElements.BarDataFigure

from ...common.common_figure_materials import CommonColorDict as CCD

TISSUE_SHORT = {
    'cortex': 'Ctx', 'hippocampus': 'Hipp', 'hypothalamus': 'Hyp',
    'liver': 'Liv', 'kidney': 'Kid', 'fat': 'Fat',
    'gastrocnemius': 'Gas', 'quadriceps': 'Quad', 'ileum': 'Ile',
    'serum': 'Ser',
}


class CrossTissueBarDataFigure(BarDataFigure):
    """跨组织分组柱状图的数据图。

    figure_data_parameter_dict:
        metabolites: list[str]  — 代谢物 ID 列表
        values: dict  — {met_id: {tissue: {group: array}}}
        metabolite_names: dict  — {met_id: display_name}
        tissue_order: list[str]  — 组织排列顺序
        group_colors: dict  — {group: color_hex}
        group_names: tuple  — (group1, group2)
    """

    def __init__(self, figure_data_parameter_dict, bottom_left: Vector, size: Vector, **kwargs):
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

        g1_name, g2_name = group_names
        n_met = len(metabolites)
        if n_met == 0:
            bar_param = {
                ParameterName.ax_bottom_left_list: [Vector(0, 0)],
                ParameterName.ax_size_list: [Vector(1, 1)],
                ParameterName.color_dict: None,
                ParameterName.data_nested_list: [({}, {})],
                ParameterName.array_len_list: [1],
                ParameterName.figure_config_dict: {},
                ParameterName.x_tick_labels_list: [['']],
                ParameterName.y_label_list: [''],
                ParameterName.legend: False,
            }
            super().__init__(bar_param, bottom_left, size, **kwargs)
            return

        color_dict = {
            g1_name: group_colors.get(g1_name, CCD.gray),
            g2_name: group_colors.get(g2_name, CCD.teal),
        }

        data_nested = []
        x_tick_labels_list = []
        array_len_list = []
        y_label_list = []
        subplot_name_list = []
        ax_bottom_left_list = []
        ax_size_list = []

        ncols = min(4, max(n_met, 1))
        nrows = (n_met + ncols - 1) // ncols
        gap_x = 0.015
        gap_y = 0.03
        ax_w = (1.0 - (ncols - 1) * gap_x) / ncols
        ax_h = (1.0 - (nrows - 1) * gap_y) / nrows

        for idx, met_id in enumerate(metabolites):
            met_values = values.get(met_id, {})
            available_tissues = [t for t in tissue_order if t in met_values]
            n_tissues = len(available_tissues)
            if n_tissues == 0:
                continue

            g1_means = np.array([
                np.mean(met_values[t][g1_name]) for t in available_tissues
                if g1_name in met_values.get(t, {})])
            g2_means = np.array([
                np.mean(met_values[t][g2_name]) for t in available_tissues
                if g2_name in met_values.get(t, {})])
            g1_stds = np.array([
                np.std(met_values[t][g1_name], ddof=1) if len(met_values[t][g1_name]) > 1 else 0
                for t in available_tissues
                if g1_name in met_values.get(t, {})])
            g2_stds = np.array([
                np.std(met_values[t][g2_name], ddof=1) if len(met_values[t][g2_name]) > 1 else 0
                for t in available_tissues
                if g2_name in met_values.get(t, {})])

            mid_dict = {g1_name: g1_means, g2_name: g2_means}
            error_dict = {g1_name: g1_stds, g2_name: g2_stds}
            data_nested.append((mid_dict, error_dict))
            x_tick_labels_list.append([TISSUE_SHORT.get(t, t[:4]) for t in available_tissues])
            array_len_list.append(n_tissues)

            row, col = divmod(idx, ncols)
            ax_bottom_left_list.append(
                Vector(col * (ax_w + gap_x), (nrows - 1 - row) * (ax_h + gap_y)))
            ax_size_list.append(Vector(ax_w, ax_h))
            y_label_list.append('Rel. Abund.' if col == 0 else None)
            subplot_name_list.append(metabolite_names.get(met_id, str(met_id)))

        bar_param = {
            ParameterName.ax_bottom_left_list: ax_bottom_left_list,
            ParameterName.ax_size_list: ax_size_list,
            ParameterName.color_dict: color_dict,
            ParameterName.data_nested_list: data_nested,
            ParameterName.array_len_list: array_len_list,
            ParameterName.figure_config_dict: {
                ParameterName.column_width: 0.6,
                ParameterName.edge: 0.1,
                ParameterName.error_bar_param_dict: DataFigureConfig.common_error_bar_param_dict_generator(),
                ParameterName.subplot_name_text_format_dict: {
                    **DataFigureConfig.common_subplot_text_format_dict_generator(),
                    ParameterName.font_size: 7,
                },
            },

            ParameterName.subplot_name_list: subplot_name_list,
            ParameterName.text_axis_loc_pair: Vector(0.5, 0.9),
            ParameterName.x_tick_labels_list: x_tick_labels_list,
            ParameterName.y_label_list: y_label_list,
            ParameterName.legend: False,
        }
        super().__init__(bar_param, bottom_left, size, **kwargs)


class CrossTissueBarFigure(CompositeFigure):
    """跨组织柱状图容器。

    figure_data_parameter_dict:
        同 CrossTissueBarDataFigure
    """

    subplot_width = 0.22
    subplot_height = 0.15
    subplot_gap_x = 0.03
    subplot_gap_y = 0.04

    def __init__(self, figure_data_parameter_dict, **kwargs):
        metabolites = figure_data_parameter_dict.get('metabolites', [])
        n_met = len(metabolites)
        ncols = min(4, max(n_met, 1))
        nrows = (n_met + ncols - 1) // ncols
        total_width = ncols * self.subplot_width + (ncols - 1) * self.subplot_gap_x
        total_height = nrows * self.subplot_height + (nrows - 1) * self.subplot_gap_y

        bar_obj = CrossTissueBarDataFigure(
            figure_data_parameter_dict=dict(figure_data_parameter_dict),
            bottom_left=Vector(0, 0),
            size=Vector(total_width, total_height),
        )

        subfigure_element_dict = {'bar': {'bar': bar_obj}}
        self.total_width = total_width
        self.height_to_width_ratio = total_height / total_width
        super().__init__(
            subfigure_element_dict,
            Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)
