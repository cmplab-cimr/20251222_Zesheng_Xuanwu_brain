from .config import Vector, ParameterName as GeneralParameterName, Elements as GeneralElements, TransformDict


class FigurePlotting(object):
    def __init__(self, ParameterName, Elements):
        self.ParameterName = ParameterName
        self.Figure = Elements.Figure
        self.Elements = Elements
        self.default_figure_size = Vector(8.5, 8.5)

    def _common_draw_function(
            self, figure_name, target_obj, figure_output_direct, figure_size=None):
        if figure_size is None:
            figure_size = self.default_figure_size
        if isinstance(target_obj, (tuple, list)):
            other_obj_list = target_obj
        else:
            other_obj_list = [target_obj]
        current_figure = self.Figure(
            figure_name, other_obj_list=other_obj_list, figure_size=figure_size,
            figure_output_direct=figure_output_direct)
        current_figure.draw()

    def violin_box_distribution_function(
            self, figure_output_direct, figure_name,
            complete_data_dict, figure_type):
        pass

    def metabolic_flux_model_function(
            self, figure_output_direct, figure_name,
            input_metabolite_set, c13_labeling_metabolite_set, mid_data_metabolite_set, mixed_mid_data_metabolite_set,
            biomass_metabolite_set, boundary_flux_set, current_reaction_value_dict=None, infusion=False,
            figure_size=None):
        ParameterName = self.ParameterName
        experimental_setting_dict = {
            ParameterName.input_metabolite_set: input_metabolite_set,
            ParameterName.c13_labeling_metabolite_set: c13_labeling_metabolite_set,
            ParameterName.mid_data_metabolite_set: mid_data_metabolite_set,
            ParameterName.mixed_mid_data_metabolite_set: mixed_mid_data_metabolite_set,
            ParameterName.biomass_metabolite_set: biomass_metabolite_set,
            ParameterName.boundary_flux_set: boundary_flux_set,
        }
        metabolic_network_config_dict = {
            ParameterName.bottom_left_offset: Vector(0.1, 0.1),
            ParameterName.scale: 0.8,
            ParameterName.infusion: infusion,
            **experimental_setting_dict
        }
        if current_reaction_value_dict is not None:
            metabolic_network_config_dict.update({
                ParameterName.reaction_raw_value_dict: current_reaction_value_dict,
                ParameterName.visualize_flux_value: ParameterName.transparency
            })
        metabolic_network_obj = self.Elements.LianfengMetabolicNetwork(**metabolic_network_config_dict)
        self._common_draw_function(
            figure_name, metabolic_network_obj, figure_output_direct, figure_size)

    def mid_prediction_function(
            self, data_name, result_label, mid_name_list, output_direct, figure_config_dict, figure_size=None):
        ParameterName = self.ParameterName
        current_mid_comparison_figure_config_dict = {
            ParameterName.bottom_left_offset: Vector(0.15, 0.05),
            ParameterName.scale: 0.7,
            ParameterName.figure_data_parameter_dict: {
                ParameterName.legend: False,
                ParameterName.data_name: data_name,
                ParameterName.result_label: result_label,
                ParameterName.mid_name_list: mid_name_list,
            },
        }
        try:
            new_figure_data_parameter_dict = figure_config_dict.pop(ParameterName.figure_data_parameter_dict)
        except KeyError:
            new_figure_data_parameter_dict = {}
        current_mid_comparison_figure_config_dict[ParameterName.figure_data_parameter_dict].update(
            new_figure_data_parameter_dict)
        current_mid_comparison_figure_config_dict.update(figure_config_dict)
        mid_comparison_figure = self.Elements.MIDComparisonGridBarWithLegendDataFigure(
            **current_mid_comparison_figure_config_dict)
        self._common_draw_function(
            result_label, mid_comparison_figure, output_direct, figure_size)

    def mid_comparison_plotting(
            self, figure_output_direct, figure_name,
            complete_data_dict, target_emu_name_nested_list, target_row_num, target_col_num,
            error_bar_data_dict=None, color_dict=None, title_dict=None,
            y_lim=(0, 1), y_label=None, x_label_list=None, figure_size=None, legend=False, legend_color_dict=None,
            supplementary_text_list=None, supplementary_text_loc_list=None, **kwargs):
        ParameterName = self.ParameterName
        if title_dict is None:
            title_dict = TransformDict()
        if target_row_num > 6:
            low_height = True
        else:
            low_height = False
        if supplementary_text_loc_list is not None:
            supplementary_text_loc_vector_list = [
                Vector(*text_loc) for text_loc in supplementary_text_loc_list]
        else:
            supplementary_text_loc_vector_list = None

        mid_comparison_grid_figure_parameter_dict = {
            ParameterName.figure_data: (complete_data_dict, error_bar_data_dict),
            ParameterName.color_dict: color_dict,
            ParameterName.common_y_lim: y_lim,
            ParameterName.legend: legend,
            ParameterName.low_height: low_height,

            ParameterName.mid_name_list: target_emu_name_nested_list,
            ParameterName.x_tick_labels_list: x_label_list,
            ParameterName.subplot_name_dict: title_dict,
            ParameterName.supplementary_text_list: supplementary_text_list,
            ParameterName.supplementary_text_loc_list: supplementary_text_loc_vector_list,

            ParameterName.x_tick_label_format_dict: {
                ParameterName.font_size: 5,
            },
            ParameterName.subplot_name_text_format_dict: {
                ParameterName.font_size: 8,
            },
            ParameterName.supplementary_text_format_dict: {
                ParameterName.font_size: 5,
                ParameterName.width: 0.5,
                ParameterName.height: 0.5,
            }
        }
        if y_label is not None:
            mid_comparison_grid_figure_parameter_dict[ParameterName.common_y_label] = y_label
        if legend_color_dict is not None:
            mid_comparison_grid_figure_parameter_dict[ParameterName.legend_color_dict] = legend_color_dict
        mid_comparison_grid_figure_config_dict = {
            ParameterName.bottom_left_offset: Vector(0.1, 0.05),
            ParameterName.scale: 0.8,
            ParameterName.figure_data_parameter_dict: mid_comparison_grid_figure_parameter_dict,
        }

        mid_comparison_obj = self.Elements.MIDComparisonGridBarWithLegendDataFigure(
            **mid_comparison_grid_figure_config_dict)
        self._common_draw_function(
            figure_name, mid_comparison_obj, figure_output_direct, figure_size)

    def multi_tumor_figure_plotting(
            self, data_name, flux_location_nested_list, output_direct, figure_size):
        ParameterName = self.ParameterName

        figure_data_parameter_dict = {
            ParameterName.data_name: ParameterName.multiple_tumor,
            ParameterName.comparison_name: '',
            ParameterName.flux_name_list: flux_location_nested_list,
            ParameterName.color_dict: {
                'kidney': 'blue',
                'lung': 'orange',
                'brain': 'purple',
            }
        }
        loss_grid_comparison_figure = self.Elements.FluxComparisonScatterWithTitle(**{
            ParameterName.total_width: 0.4,
            ParameterName.bottom_left_offset: Vector(0.15, 0.15),
            ParameterName.scale: 2,
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
        })
        figure_name = f'grid_flux_comparison_{data_name}'
        self._common_draw_function(
            figure_name, loss_grid_comparison_figure, output_direct, figure_size)

    def distribution_comparison_plotting(
            self, figure_output_direct, figure_name,
            complete_data_dict, flux_name_nested_list, display_flux_dict,
            name_dict=None, color_dict=None, figure_size=None, legend=False):
        ParameterName = self.ParameterName

        figure_data_parameter_dict = {
            ParameterName.figure_data: complete_data_dict,
            ParameterName.name_dict: name_dict,
            ParameterName.legend_color_dict: color_dict,
            ParameterName.color_dict: color_dict,
            ParameterName.flux_name_list: flux_name_nested_list,
            ParameterName.display_flux_name_dict: display_flux_dict,
            ParameterName.low_height: True,
            ParameterName.legend: legend,
            ParameterName.figure_config_dict: {
                ParameterName.column_width: 0.85,
                ParameterName.error_bar_param_dict: {
                    ParameterName.cap_size: 0.5,
                    ParameterName.edge_width: 0.1,
                },
                ParameterName.bar_param_dict: {
                    ParameterName.alpha: 1,
                }
            }
        }
        flux_comparison_obj = self.Elements.FluxComparisonGridBarWithLegendDataFigure(**{
            ParameterName.bottom_left_offset: Vector(0.1, 0.05),
            ParameterName.scale: 0.85,
            ParameterName.figure_data_parameter_dict: figure_data_parameter_dict,
        })
        self._common_draw_function(
            figure_name, flux_comparison_obj, figure_output_direct, figure_size)


figure_plotting = FigurePlotting(GeneralParameterName, GeneralElements)
