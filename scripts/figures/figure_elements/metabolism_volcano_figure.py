from ..common.config import (
    ParameterName, Vector, FontWeight, CompositeFigure, TextBox, merge_axis_format_dict,
    default_parameter_extract, GeneralElements, BasicFigureData, FigureDataKeywords, ColorConfig,
    CommonElementConfig, np, DataFigureConfig, LineStyle)


ScatterDataFigure = GeneralElements.ScatterDataFigure


mid_comparison_name_dict = {
    ParameterName.optimized: 'Optimized MID',
    ParameterName.experimental: 'Target experimental MID',
}
mid_comparison_color_dict = {
    ParameterName.optimized: ColorConfig.dark_blue,
    ParameterName.experimental: ColorConfig.orange,
}

class VolcanoScatterConfig(object):
    default_color_dict = {
        ParameterName.positive: ColorConfig.normal_blue,
        ParameterName.negative: ColorConfig.orange,
        ParameterName.medium_data: ColorConfig.light_gray,
    }


class VolcanoScatterDataFigure(ScatterDataFigure):
    def __init__(
            self, figure_data_parameter_dict, bottom_left: Vector, size: Vector, **kwargs):
        ax_total_bottom_left = Vector(0, 0)
        ax_total_size = Vector(1, 1) - ax_total_bottom_left

        (
            log_fold_data_array, log_p_value_data_array, log_fold_threshold, log_p_value_threshold
        ) = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_data, '', pop=True, force=True)
        color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.color_dict, VolcanoScatterConfig.default_color_dict, pop=True)
        x_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.common_x_label, 'log2(Fold)', pop=True)
        y_label = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.common_y_label, 'log10(Pvalue)', pop=True)
        data_name_list = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.flux_name_list, None, pop=True)
        if data_name_list is not None:
            data_name_array = np.array(data_name_list)
        else:
            data_name_array = None

        marker_size = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.marker_size, 8, pop=True)
        new_figure_config_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_config_dict, {}, pop=True)

        supplementary_text_format_dict = {
            **DataFigureConfig.common_supplementary_text_config_dict,
            ParameterName.font_size: 5,
            ParameterName.width: 0.1,
            ParameterName.height: 0.05,
            **default_parameter_extract(
                figure_data_parameter_dict, ParameterName.supplementary_text_format_dict,
                {}, pop=True)
        }

        flux_title_text_format_dict = {
            **DataFigureConfig.common_subplot_text_format_dict_generator(),
            ParameterName.font_size: 7,
        }

        figure_config_dict = {
            ParameterName.y_label_format_dict: merge_axis_format_dict(
                {}, {
                    ParameterName.axis_label_distance: 0.04
                }, new_figure_config_dict, ParameterName.y_label_format_dict),
            ParameterName.x_tick_label_format_dict: merge_axis_format_dict(
                {}, {
                    ParameterName.axis_tick_label_distance: 0.008
                }, new_figure_config_dict, ParameterName.x_tick_label_format_dict),
            ParameterName.x_label_format_dict: merge_axis_format_dict(
                {}, {
                    ParameterName.axis_label_distance: 0.025
                }, new_figure_config_dict, ParameterName.x_label_format_dict),
            ParameterName.subplot_name_text_format_dict: flux_title_text_format_dict,
            ParameterName.line_param_dict: {
                **DataFigureConfig.common_line_param_dict_generator(),
                ParameterName.z_order: DataFigureConfig.line_z_order,
                ParameterName.edge_style: LineStyle.dash,
            },
            ParameterName.supplementary_text_format_dict: supplementary_text_format_dict
        }

        x_value_array_list = []
        y_value_array_list = []
        marker_color_list = []
        supplementary_text_list = []
        supplementary_text_loc_list = []

        x_min = np.min(log_fold_data_array)
        x_max = np.max(log_fold_data_array)
        y_max = np.max(log_p_value_data_array)
        x_lim = Vector(x_min, x_max) + Vector(-1, 1) * 0.1 * (x_max - x_min)
        y_lim = (-0.9, 1.1 * y_max)
        positive_fold_data_bool_array = (
            (log_fold_data_array > log_fold_threshold) &
            (log_p_value_data_array > log_p_value_threshold))
        negative_fold_data_bool_array = (
            (log_fold_data_array < -log_fold_threshold) &
            (log_p_value_data_array > log_p_value_threshold))
        insignificant_data_bool_array = (~positive_fold_data_bool_array) & (~negative_fold_data_bool_array)
        for data_label, current_bool_array in zip(
                [ParameterName.positive, ParameterName.negative, ParameterName.medium_data,],
                [positive_fold_data_bool_array, negative_fold_data_bool_array, insignificant_data_bool_array]):
            current_color = color_dict[data_label]
            current_fold_data_array = log_fold_data_array[current_bool_array]
            current_p_value_array = log_p_value_data_array[current_bool_array]
            x_value_array_list.append(current_fold_data_array)
            y_value_array_list.append(current_p_value_array)
            array_num = np.count_nonzero(current_bool_array)
            marker_color_list.extend([current_color] * array_num)
            if data_name_array is not None and (
                    data_label == ParameterName.positive or data_label == ParameterName.negative):
                supplementary_text_list.extend(data_name_array[current_bool_array])
                current_data_loc_array = np.vstack(
                    [current_fold_data_array, current_p_value_array]).T
                current_axis_loc_array = (
                    current_data_loc_array - np.array([[x_lim[0], y_lim[0]]])
                ) / np.array([[(x_lim[1] - x_lim[0]), (y_lim[1] - y_lim[0])]])
                supplementary_text_loc_list.extend(current_axis_loc_array + np.array([[0, 0.02]]))

        current_data_dict = {
            ParameterName.x_value_array: np.concatenate(x_value_array_list),
            ParameterName.y_value_array: (np.concatenate(y_value_array_list), None),
            ParameterName.marker_size: marker_size,
            ParameterName.marker_color: marker_color_list,  # To avoid warning from matplotlib
            ParameterName.scatter_param_dict: {
                ParameterName.alpha: 0.8,
                ParameterName.z_order: DataFigureConfig.normal_figure_element_z_order
            },
        }
        common_dash_line_config_dict = {
            **DataFigureConfig.common_line_param_dict_generator(),
            ParameterName.edge_color: ColorConfig.gray,
            ParameterName.z_order: DataFigureConfig.line_z_order,
            ParameterName.edge_style: LineStyle.dash,
        }

        scatter_line_pair_list = [[
            [x_lim[0], x_lim[1],],
            [log_p_value_threshold, log_p_value_threshold,],
            dict(common_dash_line_config_dict),
        ],[
            [log_fold_threshold, log_fold_threshold,],
            [y_lim[0], y_lim[1],],
            dict(common_dash_line_config_dict),
        ],[
            [-log_fold_threshold, -log_fold_threshold,],
            [y_lim[0], y_lim[1],],
            dict(common_dash_line_config_dict),
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
            **figure_data_parameter_dict
        }

        super().__init__(figure_data_parameter_dict, bottom_left, size, **kwargs)


class MetabolismVolcanoFigure(CompositeFigure):
    height_to_width_ratio = 0.8
    legend_height = 0.06
    title_height = 0.05

    def __init__(
            self, figure_data_parameter_dict, **kwargs):
        total_width = self.total_width
        figure_height = 0.8
        figure_title = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_title, '', pop=True)

        figure_bottom_line = 0.05
        figure_left = 0.1
        figure_width = total_width - figure_left
        title_height = self.title_height
        title_center_x = figure_left + figure_width / 2
        title_center_y = figure_height + title_height / 2
        total_height = figure_height + self.title_height
        self.total_height = total_height
        self.height_to_width_ratio = total_height / total_width

        volcano_config_dict = {
            ParameterName.bottom_left: (figure_left, figure_bottom_line),
            ParameterName.size: Vector(figure_width, figure_height),
            ParameterName.figure_data_parameter_dict: {
                **figure_data_parameter_dict,
            },
        }

        main_title_text_config_dict = {
            **CommonElementConfig.common_text_config,
            ParameterName.font_weight: FontWeight.bold,
            ParameterName.text_box: False,
            ParameterName.string: figure_title,
            ParameterName.font_size: 17,
            ParameterName.width: total_width,
            ParameterName.height: title_height,
            ParameterName.center: Vector(title_center_x, title_center_y)
        }

        subfigure_element_dict = {
            'mid_comparison': {
                'mid_comparison': VolcanoScatterDataFigure(**volcano_config_dict)},
            'title': {
                'title': TextBox(**main_title_text_config_dict)},
        }
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height),
            background=False, **kwargs)




