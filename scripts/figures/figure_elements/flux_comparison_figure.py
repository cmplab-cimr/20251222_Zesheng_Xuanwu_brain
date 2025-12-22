from ..common.config import DataFigureConfig, ParameterName, Vector, FontWeight, CompositeFigure, \
    default_parameter_extract
from ..common.common_figure_materials import CommonColorDict, CommonFigureString
from .violin_box_data_figure import CommonComparisonViolinBoxDataFigure, \
    LowRowHeightCommonComparisonViolinBoxDataFigure, ComparisonType


class FluxComparisonGridBarWithLegendDataFigure(CompositeFigure):
    height_to_width_ratio = 0.8
    legend_height = 0.06

    @staticmethod
    def calculate_height(self, scale, DataFigureClass, total_row_num, legend=False):
        figure_height = DataFigureClass.calculate_height(DataFigureClass, total_row_num)
        if legend:
            total_height = figure_height + self.legend_height
        else:
            total_height = figure_height
        return figure_height, total_height

    def __init__(
            self, figure_data_parameter_dict, total_width=1, scale=1, **kwargs):
        self.total_width = total_width
        flux_name_nested_list = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.flux_name_list, None, force=True, pop=True)
        legend = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.legend, False, force=True, pop=True)
        legend_name_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.name_dict, None, pop=True, force=True)
        color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.color_dict, None, pop=True, force=True)
        legend_color_dict = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.legend_color_dict, color_dict)
        total_row_num = len(flux_name_nested_list)
        low_height = default_parameter_extract(
            figure_data_parameter_dict, ParameterName.low_height, False, pop=True)
        if low_height:
            DataFigureClass = LowRowHeightCommonComparisonViolinBoxDataFigure
        else:
            DataFigureClass = CommonComparisonViolinBoxDataFigure
        figure_height, total_height = FluxComparisonGridBarWithLegendDataFigure.calculate_height(
            self, scale, DataFigureClass, total_row_num, legend)
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

        flux_comparison_config_dict = {
            ParameterName.bottom_left: (0, bottom_line),
            ParameterName.size: Vector(total_width, figure_height),
            ParameterName.figure_data_parameter_dict: {
                ParameterName.color_dict: color_dict,
                ParameterName.flux_name_list: flux_name_nested_list,
                ParameterName.legend: legend,
                ParameterName.legend_config_dict: legend_config_dict,
                ParameterName.size: Vector(total_width, figure_height),
                **figure_data_parameter_dict,
            },
        }

        subfigure_element_dict = {
            'flux_comparison': {
                'flux_comparison': DataFigureClass(
                    **flux_comparison_config_dict, comparison_type=ComparisonType.flux_comparison)},
        }
        super().__init__(
            subfigure_element_dict, Vector(0, 0), Vector(total_width, total_height), scale=scale,
            background=False, **kwargs)

