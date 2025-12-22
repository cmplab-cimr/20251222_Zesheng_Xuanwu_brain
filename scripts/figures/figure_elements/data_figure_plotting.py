from ..figure_plotting_package.common.classes import Vector, np
from ..figure_plotting_package.common.figure_data_format import check_and_mkdir_of_direct


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
        check_and_mkdir_of_direct(figure_output_direct)
        current_figure = self.Figure(
            figure_name, other_obj_list=other_obj_list, figure_size=figure_size,
            figure_output_direct=figure_output_direct)
        current_figure.draw()

    def _base_log(self, data_vector, base):
        if base == 2:
            result = np.log2(data_vector)
        elif base == 10:
            result = np.log10(data_vector)
        else:
            raise ValueError()
        if np.any(np.isnan(result)):
            raise ValueError()
        return result

    def volcano_plot(
            self, figure_output_direct, figure_name, figure_title,
            fold_data_array, p_value_data_array, data_name_list, fold_log_threshold, p_value_log_threshold,
            fold_data_log_base_num=2, p_value_data_log_base_num=10, figure_size=None):
        ParameterName = self.ParameterName

        log_fold_data_array = self._base_log(fold_data_array, fold_data_log_base_num)
        log_p_value_data_array = -self._base_log(p_value_data_array, fold_data_log_base_num)
        figure_data_tuple = (
            log_fold_data_array, log_p_value_data_array, fold_log_threshold, p_value_log_threshold)

        metabolic_network_config_dict = {
            ParameterName.bottom_left_offset: Vector(0.1, 0.1),
            ParameterName.scale: 0.8,
            ParameterName.figure_data_parameter_dict: {
                ParameterName.figure_title: figure_title,
                ParameterName.figure_data: figure_data_tuple,
                ParameterName.common_x_label: f'log{fold_data_log_base_num}(Fold)',
                ParameterName.common_y_label: f'log{p_value_data_log_base_num}(Pvalue)',
                ParameterName.flux_name_list: data_name_list,
            }
        }
        volcano_plotting_obj = self.Elements.MetabolismVolcanoFigure(**metabolic_network_config_dict)
        self._common_draw_function(
            figure_name, volcano_plotting_obj, figure_output_direct, figure_size)

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

