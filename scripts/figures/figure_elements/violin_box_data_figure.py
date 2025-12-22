from ..common.config import DataFigureConfig, ParameterName, Vector, Keywords, \
    np, default_parameter_extract, VerticalAlignment, merge_axis_format_dict, ColorConfig, \
    generate_violin_config_dict, merge_complete_config_dict, GeneralElements

ViolinBoxDataFigure = GeneralElements.ViolinBoxDataFigure


class ComparisonType(object):
    flux_comparison = 'flux_comparison'
    loss_comparison = 'loss_comparison'


def flux_comparison_generator(
        figure_data_parameter_dict, ax_total_bottom_left, ax_total_size, ax_interval,
        ax_bottom_left_list, ax_size_list, data_nested_list, positions_list,
        x_lim_list, x_label_list, x_ticks_list, x_tick_labels_list,
        y_lim_list, y_label_list, y_ticks_list, y_tick_labels_list,
        preset_x_lim_list, preset_y_lim_list, preset_y_ticks_list, common_x_label, common_y_label,
        color_list, flux_name_list, column_width, gap_inside_column, color_dict,
        display_group_name_dict, display_flux_name_dict, ):
    final_flux_comparison_data_dict = default_parameter_extract(
        figure_data_parameter_dict, ParameterName.figure_data, None, force=True, pop=True)

    flux_name_nested_list = figure_data_parameter_dict[ParameterName.flux_name_list]
    row_num = len(flux_name_nested_list)
    ax_row_size = (ax_total_size.y - (row_num - 1) * ax_interval.y) / row_num
    min_column_width = np.inf

    for row_index, row_list in enumerate(flux_name_nested_list):
        this_row_axis_num = len(row_list)
        this_row_axis_col_width = (ax_total_size.x - (this_row_axis_num - 1) * ax_interval.x) / this_row_axis_num
        this_row_bottom_left = (
            ax_total_bottom_left + Vector(0, (ax_row_size + ax_interval.y) * (row_num - row_index - 1)))
        for _ in range(this_row_axis_num):
            ax_size_list.append(Vector(this_row_axis_col_width, ax_row_size))
            ax_bottom_left_list.append(this_row_bottom_left)
            this_row_bottom_left = this_row_bottom_left + Vector(this_row_axis_col_width + ax_interval.x, 0)
        for col_index, flux_name in enumerate(row_list):
            current_ax_data_list = []
            current_ax_position_list = []
            current_ax_color_list = []
            current_flux_data_dict = final_flux_comparison_data_dict[flux_name]
            current_x_tick_labels = []
            group_class_dict = {}
            for group_data_dict in current_flux_data_dict.values():
                group_class_dict.update(group_data_dict)
            max_class_num = len(group_class_dict)
            total_group_num = len(current_flux_data_dict)
            each_column_width_inside_group = \
                column_width / (max_class_num + gap_inside_column * (max_class_num - 1))
            absolute_gap = each_column_width_inside_group * gap_inside_column
            min_column_width = np.minimum(min_column_width, each_column_width_inside_group)
            # In this part, group is different patient/sample, while class is different condition.
            for group_index, (group_name, group_data_dict) in enumerate(current_flux_data_dict.items()):
                group_location = group_index + 1
                for class_index, class_name in enumerate(group_class_dict.keys()):
                    data_array = group_data_dict[class_name]
                    current_ax_data_list.append(data_array)
                    current_ax_position_list.append(
                        group_location - column_width / 2 +
                        (each_column_width_inside_group + absolute_gap) * class_index +
                        each_column_width_inside_group / 2)
                    current_ax_color_list.append(color_dict[class_name])
                if group_name in display_group_name_dict:
                    display_group_name = display_group_name_dict[group_name]
                else:
                    display_group_name = group_name
                current_x_tick_labels.append(display_group_name)

            if flux_name in display_flux_name_dict:
                display_flux_name = display_flux_name_dict[flux_name]
            else:
                display_flux_name = flux_name
            flux_name_list.append(display_flux_name)
            x_ticks_list.append(list(range(1, total_group_num + 1)))
            data_nested_list.append(current_ax_data_list)
            positions_list.append(current_ax_position_list)
            y_tick_labels_list.append(Keywords.default)
            if preset_x_lim_list is not None:
                x_lim = preset_x_lim_list[row_index][col_index]
            else:
                x_lim = (0.5, total_group_num + 0.5)
            x_lim_list.append(x_lim)
            if preset_y_lim_list is not None:
                y_lim = preset_y_lim_list[row_index][col_index]
            else:
                y_lim = None
            y_lim_list.append(y_lim)
            if len(color_list) == 0:
                color_list.extend(current_ax_color_list)
            if row_index == row_num - 1:
                x_tick_labels_list.append(current_x_tick_labels)
                x_label_list.append(common_x_label)
            else:
                x_tick_labels_list.append(None)
                x_label_list.append(None)
            if col_index == 0:
                y_label_list.append(common_y_label)
            else:
                y_label_list.append(None)
            if preset_y_ticks_list is not None:
                y_ticks = preset_y_ticks_list[row_index][col_index]
            else:
                y_ticks = None
            y_ticks_list.append(y_ticks)
    return min_column_width


def loss_comparison_generator(
        figure_data_parameter_dict, ax_total_bottom_left, ax_total_size, ax_interval,
        ax_bottom_left_list, ax_size_list, data_nested_list, positions_list,
        x_lim_list, x_label_list, x_ticks_list, x_tick_labels_list,
        y_lim_list, y_label_list, y_ticks_list, y_tick_labels_list,
        preset_x_lim_list, preset_y_lim_list, preset_y_ticks_list, preset_y_tick_labels_list,
        common_x_label, common_y_label,
        color_list, column_width, gap_inside_column, color_dict,
        display_group_name_dict, ):
    loss_data_dict = default_parameter_extract(
        figure_data_parameter_dict, ParameterName.figure_data, None, force=True, pop=True)
    result_label_layout_list = default_parameter_extract(
        figure_data_parameter_dict, ParameterName.result_label_layout_list, None)
    if result_label_layout_list is None:
        result_label_layout_list = [['']]
        loss_data_dict = {(0, 0): loss_data_dict}
        row_num = 1
        if preset_y_lim_list is not None:
            preset_y_lim_list = [[preset_y_lim_list]]
        if preset_y_ticks_list is not None:
            preset_y_ticks_list = [[preset_y_ticks_list]]
        if preset_y_tick_labels_list is not Keywords.default:
            preset_y_tick_labels_list = [[preset_y_tick_labels_list]]
    else:
        row_num = len(result_label_layout_list)

    ax_row_size = (ax_total_size.y - (row_num - 1) * ax_interval.y) / row_num
    min_column_width = np.inf
    ax_total_width = 1

    for row_index, row_list in enumerate(result_label_layout_list):
        this_row_axis_num = len(row_list)
        this_row_axis_col_width = (ax_total_width - (this_row_axis_num - 1) * ax_interval.x) / this_row_axis_num
        this_row_bottom_left = ax_total_bottom_left + \
                               Vector(0, (ax_row_size + ax_interval.y) * (row_num - row_index - 1))
        for _ in range(this_row_axis_num):
            ax_size_list.append(Vector(this_row_axis_col_width, ax_row_size))
            ax_bottom_left_list.append(this_row_bottom_left)
            this_row_bottom_left = this_row_bottom_left + Vector(this_row_axis_col_width + ax_interval.x, 0)
        for col_index, result_label in enumerate(row_list):
            current_ax_data_list = []
            current_ax_position_list = []
            current_ax_color_list = []
            current_loss_data_dict = loss_data_dict[(row_index, col_index)]
            current_x_tick_labels = []
            group_class_dict = {}
            for group_data_dict in current_loss_data_dict.values():
                group_class_dict.update(group_data_dict)
            max_class_num = len(group_class_dict)
            total_group_num = len(current_loss_data_dict)
            each_column_width_inside_group = \
                column_width / (max_class_num + gap_inside_column * (max_class_num - 1))
            absolute_gap = each_column_width_inside_group * gap_inside_column
            min_column_width = np.minimum(min_column_width, each_column_width_inside_group)
            # In this part, group is different patient/sample, while class is different condition.
            for group_index, (group_name, group_data_dict) in enumerate(current_loss_data_dict.items()):
                group_location = group_index + 1
                for class_index, class_name in enumerate(group_class_dict.keys()):
                    data_array = group_data_dict[class_name]
                    current_ax_data_list.append(data_array)
                    current_ax_position_list.append(
                        group_location - column_width / 2 +
                        (each_column_width_inside_group + absolute_gap) * class_index +
                        each_column_width_inside_group / 2)
                    current_ax_color_list.append(color_dict[class_name])
                if group_name in display_group_name_dict:
                    display_group_name = display_group_name_dict[group_name]
                else:
                    display_group_name = group_name
                current_x_tick_labels.append(display_group_name)

            x_ticks_list.append(list(range(1, total_group_num + 1)))
            data_nested_list.append(current_ax_data_list)
            positions_list.append(current_ax_position_list)
            if preset_y_tick_labels_list != Keywords.default:
                y_tick_labels_list.append(preset_y_tick_labels_list[row_index][col_index])
            else:
                y_tick_labels_list.append(Keywords.default)
            if preset_x_lim_list is not None:
                x_lim = preset_x_lim_list[row_index][col_index]
            else:
                x_lim = (0.5, total_group_num + 0.5)
            x_lim_list.append(x_lim)
            if preset_y_lim_list is not None:
                y_lim = preset_y_lim_list[row_index][col_index]
            else:
                y_lim = None
            y_lim_list.append(y_lim)
            if len(color_list) == 0:
                color_list.extend(current_ax_color_list)
            if row_index == row_num - 1:
                x_tick_labels_list.append(current_x_tick_labels)
                x_label_list.append(common_x_label)
            else:
                x_tick_labels_list.append(None)
                x_label_list.append(None)
            if col_index == 0:
                y_label_list.append(common_y_label)
            else:
                y_label_list.append(None)
            if preset_y_ticks_list is not None:
                y_ticks = preset_y_ticks_list[row_index][col_index]
            else:
                y_ticks = None
            y_ticks_list.append(y_ticks)
    return min_column_width


class CommonComparisonViolinBoxDataFigure(ViolinBoxDataFigure):
    ax_interval = Vector(0.015, 0.03)        # (horizontal, vertical)
    each_row_figure_height = 0.4

    @staticmethod
    def calculate_height(self, row_num, each_row_figure_height=None):
        if each_row_figure_height is None:
            each_row_figure_height = self.each_row_figure_height
        return each_row_figure_height * row_num + (row_num - 1) * self.ax_interval.y

    def __init__(
            self, figure_data_parameter_dict, bottom_left: Vector, size: Vector,
            comparison_type=ComparisonType.flux_comparison, **kwargs):
        ax_total_bottom_left = Vector(0, 0)
        ax_total_size = Vector(1, 1) - ax_total_bottom_left
        ax_interval = Vector(0.045, 0.035)

        (
            common_x_label, common_y_label,
            preset_y_lim_list, preset_y_ticks_list, preset_y_tick_labels_list, preset_x_lim_list,
            display_flux_name_dict, display_group_name_dict, figure_type
        ) = default_parameter_extract(
            figure_data_parameter_dict, [
                ParameterName.common_x_label, ParameterName.common_y_label,
                ParameterName.y_lim_list, ParameterName.y_ticks_list, ParameterName.y_tick_labels_list,
                ParameterName.x_lim_list,
                ParameterName.display_flux_name_dict, ParameterName.display_group_name_dict, ParameterName.figure_type
            ], [
                None, 'Flux value',
                None, None, Keywords.default,
                None,
                {}, {}, ParameterName.box
            ], pop=True)
        color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.color_dict, None, force=True, pop=True)
        each_row_figure_height = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.each_row_figure_height, self.each_row_figure_height)
        self.each_row_figure_height = each_row_figure_height

        text_axis_loc_pair = Vector(0.5, 1.07)
        common_line_width = 0.5
        column_width = 0.7
        gap_inside_column = 0.1
        box_body_alpha = 0.3
        flux_title_text_format_dict = {
            **DataFigureConfig.common_subplot_text_format_dict_generator(),
            ParameterName.font_size: 9,
        }

        data_nested_list = []
        positions_list = []
        x_label_list = []
        x_ticks_list = []
        x_tick_labels_list = []
        y_label_list = []
        y_ticks_list = []
        y_tick_labels_list = []
        x_lim_list = []
        y_lim_list = []
        ax_bottom_left_list = []
        ax_size_list = []
        flux_name_list = []
        color_list = []

        if comparison_type == ComparisonType.flux_comparison:
            min_column_width = flux_comparison_generator(
                figure_data_parameter_dict, ax_total_bottom_left, ax_total_size, ax_interval,
                ax_bottom_left_list, ax_size_list, data_nested_list, positions_list,
                x_lim_list, x_label_list, x_ticks_list, x_tick_labels_list,
                y_lim_list, y_label_list, y_ticks_list, y_tick_labels_list,
                preset_x_lim_list, preset_y_lim_list, preset_y_ticks_list, common_x_label, common_y_label,
                color_list, flux_name_list, column_width, gap_inside_column, color_dict,
                display_group_name_dict, display_flux_name_dict, )
        elif comparison_type == ComparisonType.loss_comparison:
            min_column_width = loss_comparison_generator(
                figure_data_parameter_dict, ax_total_bottom_left, ax_total_size, ax_interval,
                ax_bottom_left_list, ax_size_list, data_nested_list, positions_list,
                x_lim_list, x_label_list, x_ticks_list, x_tick_labels_list,
                y_lim_list, y_label_list, y_ticks_list, y_tick_labels_list,
                preset_x_lim_list, preset_y_lim_list, preset_y_ticks_list, preset_y_tick_labels_list,
                common_x_label, common_y_label,
                color_list, column_width, gap_inside_column, color_dict,
                display_group_name_dict, )
        else:
            raise ValueError()

        general_figure_config_dict = {
            ParameterName.x_tick_label_format_dict: {
                ParameterName.font_size: 7
            },
            ParameterName.x_label_format_dict: {
                ParameterName.axis_label_distance: 0.03,
                ParameterName.font_size: 10
            },
            ParameterName.y_label_format_dict: {
                ParameterName.axis_label_distance: 0.04,
                ParameterName.font_size: 10
            },
            ParameterName.box_violin_config_dict: generate_violin_config_dict(
                min_column_width, box_body_alpha, common_line_width, color_list, color_list),
            ParameterName.subplot_name_text_format_dict: flux_title_text_format_dict,
        }
        specific_figure_config_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.figure_config_dict, {}, pop=True)
        figure_config_dict = merge_complete_config_dict(general_figure_config_dict, specific_figure_config_dict)

        figure_data_parameter_dict = {
            ParameterName.figure_type: figure_type,
            ParameterName.ax_bottom_left_list: ax_bottom_left_list,
            ParameterName.ax_size_list: ax_size_list,
            ParameterName.legend_center: None,
            ParameterName.legend_area_size: None,
            ParameterName.color_dict: color_dict,
            ParameterName.name_dict: None,
            ParameterName.legend: False,
            ParameterName.figure_config_dict: figure_config_dict,

            ParameterName.data_nested_list: data_nested_list,
            ParameterName.positions_list: positions_list,
            ParameterName.cutoff: None,
            ParameterName.emphasized_flux_list: None,
            ParameterName.x_lim_list: x_lim_list,
            ParameterName.y_lim_list: y_lim_list,
            ParameterName.x_label_list: x_label_list,
            ParameterName.x_ticks_list: x_ticks_list,
            ParameterName.y_label_list: y_label_list,
            ParameterName.x_tick_labels_list: x_tick_labels_list,
            ParameterName.y_ticks_list: y_ticks_list,
            ParameterName.y_tick_labels_list: y_tick_labels_list,
            ParameterName.subplot_name_list: flux_name_list,
            ParameterName.text_axis_loc_pair: text_axis_loc_pair,
            **figure_data_parameter_dict
        }
        super().__init__(figure_data_parameter_dict, bottom_left, size, **kwargs)


class LowRowHeightCommonComparisonViolinBoxDataFigure(CommonComparisonViolinBoxDataFigure):
    each_row_figure_height = 0.1
