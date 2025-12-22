from ..common.config import DataFigureConfig, ParameterName, Vector, FontWeight, CompositeFigure, \
    default_parameter_extract, GeneralElements

MIDComparisonGridBarDataFigure = GeneralElements.MIDComparisonGridBarDataFigure


class LowHeightMIDComparisonGridBarDataFigure(MIDComparisonGridBarDataFigure):
    each_row_figure_height = 0.09


class MIDComparisonGridBarWithLegendDataFigure(CompositeFigure):
    height_to_width_ratio = 0.8
    legend_height = 0.06

    def __init__(
            self, figure_data_parameter_dict, total_width=1, scale=1, **kwargs):
        self.total_width = total_width
        mid_name_list = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.mid_name_list, None, force=True, pop=True)
        color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.color_dict, None,
            force=True, pop=True)
        legend_color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.legend_color_dict, color_dict)
        default_name_dict = {key: key for key in legend_color_dict.keys()}
        legend_name_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.name_dict, default_name_dict, pop=True)
        legend_height = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.legend_height, self.legend_height)
        self.legend_height = legend_height
        total_row_num = len(mid_name_list)
        low_height = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.low_height, False, pop=True)
        if low_height:
            figure_class = LowHeightMIDComparisonGridBarDataFigure
        else:
            figure_class = MIDComparisonGridBarDataFigure

        figure_height = figure_class.calculate_height(figure_class, total_row_num)
        total_height = figure_height + legend_height
        self.total_height = total_height
        self.height_to_width_ratio = total_height / total_width
        bottom_line = 0
        legend_bottom = figure_height + 0.005
        legend_top = total_height - 0.005
        legend_config_dict = {
            ParameterName.legend_center: Vector(0.5 * total_width, (legend_top + legend_bottom) / 2),
            ParameterName.legend_area_size: Vector(total_width, legend_top - legend_bottom),
            ParameterName.legend_color_dict: legend_color_dict,
            ParameterName.name_dict: legend_name_dict,
            ParameterName.text_config_dict: {
                ParameterName.font_size: 10,
                ParameterName.font_weight: FontWeight.bold
            }
        }

        mid_comparison_config_dict = {
            ParameterName.bottom_left: (0, bottom_line),
            ParameterName.size: Vector(total_width, figure_height),
            ParameterName.figure_data_parameter_dict: {
                ParameterName.color_dict: color_dict,
                ParameterName.mid_name_list: mid_name_list,
                ParameterName.legend: True,
                ParameterName.legend_config_dict: legend_config_dict,
                ParameterName.size: Vector(total_width, figure_height),
                **figure_data_parameter_dict,
            },
        }

        subfigure_element_dict = {
            'mid_comparison': {
                'mid_comparison': figure_class(**mid_comparison_config_dict)},
        }
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height), scale=scale,
            background=False, **kwargs)


