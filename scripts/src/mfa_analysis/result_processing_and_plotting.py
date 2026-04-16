from .config import (
    GeneralFinalResult, CoreConstants, multiple_repeat_result_processing, GeneralKeywords,
    mid_name_process, check_if_mix_flux, tissue_name_breakdown, loss_data_distribution_plotting,
    time_distribution_plotting, DataSource, FigureData, FigureDataKeywords,
    tissue_specific_name_constructor, np, reversible_flux_pair_dict_generator,
    tissue_specific_reversible_flux_title_constructor_generator, check_and_mkdir_of_direct, figure_plotting,
    Direct, experimental_data_plotting_func_template, output_raw_experimental_mid_data,
    TransformDict, output_raw_flux_data, output_predicted_mid_data, output_flux_comparison_data,
    experimental_mid_prediction_plotting_func_template, flux_comparison_plotting_func_template,
)


class FinalResult(GeneralFinalResult):
    def __init__(
            self, project_output_direct, common_data_output_direct, result_name,
            predicted_data_figure=False, suffix=None):
        super().__init__(
            project_output_direct, common_data_output_direct, result_name, suffix=suffix)
        self.mid_name_tissue_name_dict = {}
        self.flux_name_tissue_name_dict = {}
        self.predicted_data_figure = predicted_data_figure
        self.mix_prefix = f'{CoreConstants.mix_ratio_prefix}{CoreConstants.mix_flux_sep}'

    def mid_name_process(self, raw_mid_name):
        modified_mid_name = mid_name_process(raw_mid_name)
        tissue_name, pure_mid_name = tissue_name_breakdown(modified_mid_name)
        if raw_mid_name not in self.mid_name_tissue_name_dict:
            self.mid_name_tissue_name_dict[raw_mid_name] = (tissue_name, pure_mid_name)
        return modified_mid_name

    def flux_name_process(self, complete_flux_name):
        mix_flux = check_if_mix_flux(complete_flux_name)
        if mix_flux:
            tissue_name, raw_flux_name = tissue_name_breakdown(complete_flux_name[len(self.mix_prefix):])
            raw_flux_name = f'{self.mix_prefix}{raw_flux_name}'
        else:
            tissue_name, raw_flux_name = tissue_name_breakdown(complete_flux_name)
        if complete_flux_name not in self.flux_name_tissue_name_dict:
            self.flux_name_tissue_name_dict[complete_flux_name] = (tissue_name, raw_flux_name)
        return tissue_name, raw_flux_name

    def target_mid_name_process(self, raw_target_mid_name):
        new_target_mid_name = raw_target_mid_name.replace(f'{CoreConstants.mix_ratio_prefix}(', '')
        new_target_mid_name = new_target_mid_name.replace(f'{CoreConstants.convolution_id}(', '')
        new_target_mid_name = new_target_mid_name.replace(f')', '')
        modified_target_mid_name = self.mid_name_process(new_target_mid_name)
        return modified_target_mid_name

    def iteration(self, raw_result_label, process_mid_name=True, result_label_suffix_tuple=(), solver_obj=None):
        emu_sep = CoreConstants.emu_carbon_list_str_sep
        tissue_sep = CoreConstants.specific_tissue_sep

        tissue_name_replacement_pair_list = [
            ('sub_q_fat', 'white_adipose'),
            ('vastus_muscle', 'muscle'),
        ]
        (
            loss_array, solution_array, flux_name_index_dict, raw_result_information_dict, predicted_data_dict,
            raw_target_experimental_mid_data_dict, time_array
        ) = super(FinalResult, self).iteration(
            raw_result_label, process_mid_name, result_label_suffix_tuple, solver_obj=solver_obj)
        result_information_dict = {**raw_result_information_dict}
        target_experimental_mid_data_dict = {**raw_target_experimental_mid_data_dict}

        return loss_array, solution_array, flux_name_index_dict, result_information_dict, predicted_data_dict, \
            target_experimental_mid_data_dict, time_array

    def final_process(
            self, result_process_func=None, solver_dict=None, **final_process_parameters):
        result_final_process(
            self, solver_dict, **final_process_parameters)


class OutOfOrderFinalResult(FinalResult):
    def __init__(self, *args, **kwargs):
        super(OutOfOrderFinalResult, self).__init__(*args, **kwargs)

    def add_and_save_result(
            self, raw_result_label, current_result_information, result_list, flux_name_index_dict,
            target_experimental_mid_data_dict, solution_id_array=None):
        result_label = self._update_result_label(raw_result_label)
        current_solution_array, current_time_array, current_loss_array, current_predicted_dict = result_list
        self._merge_to_final_result_dict_with_solution_id(
            result_label, current_solution_array, current_time_array, current_loss_array,
            current_predicted_dict, current_result_information, flux_name_index_dict, target_experimental_mid_data_dict,
            solution_id_array)
        self._save_data(
            result_label, self.final_solution_data_dict[result_label],
            self.final_time_data_dict[result_label], self.final_loss_data_dict[result_label],
            self.final_predicted_data_dict[result_label], current_result_information,
            flux_name_index_dict, target_experimental_mid_data_dict, self.final_solution_id_array_dict[result_label])

    def load_previous_results(self, raw_result_label, mfa_tuple=None):
        return self._load_previous_results(raw_result_label, out_of_order=True, mfa_tuple=mfa_tuple)


def result_final_process(
        final_result_obj, solver_dict, result_process_information_dict_analysis=None,
        benchmark=False, average_optimized_results=None, loss_percentile=None,
        result_label_suffix_tuple=(), parallel_num=None, **other_param_dict):
    for current_result_label in solver_dict.keys():
        solver_obj = solver_dict[current_result_label]
        final_result_obj.load_current_result_label(
            current_result_label, result_label_suffix_tuple, solver_obj=solver_obj)

    final_solution_data_dict = final_result_obj.final_solution_data_dict
    final_loss_data_dict = final_result_obj.final_loss_data_dict
    final_flux_name_index_dict = final_result_obj.final_flux_name_index_dict
    final_time_data_dict = final_result_obj.final_time_data_dict
    experiment_name = final_result_obj.result_name.value
    final_predicted_data_dict = final_result_obj.final_predicted_data_dict
    final_target_experimental_mid_data_dict = final_result_obj.final_target_experimental_mid_data_dict
    final_solution_id_array_dict = final_result_obj.final_solution_id_array_dict
    mid_name_tissue_name_dict = final_result_obj.mid_name_tissue_name_dict
    final_information_dict = final_result_obj.final_information_dict

    if average_optimized_results is not None:
        (
            averaged_final_solution_data_dict, averaged_final_loss_data_dict, averaged_final_predicted_data_dict
        ) = multiple_repeat_result_processing(
            solver_dict, final_solution_data_dict, final_loss_data_dict, final_predicted_data_dict,
            repeat_division_num=average_optimized_results, loss_percentile=loss_percentile)
        loss_percentile = None
    else:
        averaged_final_solution_data_dict = final_solution_data_dict
        averaged_final_loss_data_dict = final_loss_data_dict
        averaged_final_predicted_data_dict = final_predicted_data_dict
    if not benchmark:
        final_solution_data_dict = averaged_final_solution_data_dict
        final_loss_data_dict = averaged_final_loss_data_dict
        final_predicted_data_dict = averaged_final_predicted_data_dict

    subset_index_dict = loss_data_distribution_plotting(
        experiment_name, final_loss_data_dict,
        output_direct=final_result_obj.this_result_output_direct, loss_percentile=loss_percentile)

    time_distribution_plotting(
        experiment_name, final_time_data_dict,
        output_direct=final_result_obj.this_result_output_direct, subset_index_dict=None, parallel_num=parallel_num)

    multi_tissue_experimental_mid_prediction(
        experiment_name, final_predicted_data_dict, final_target_experimental_mid_data_dict,
        solver_dict, mid_name_tissue_name_dict, final_information_dict, result_process_information_dict_analysis,
        prediction_output_direct=final_result_obj.mid_prediction_output_direct,
        subset_index_dict=subset_index_dict, **other_param_dict)

    flux_comparison_parameter_generator(
        experiment_name, averaged_final_solution_data_dict, final_flux_name_index_dict,
        final_information_dict, final_result_obj.flux_name_tissue_name_dict,
        result_process_information_dict_analysis,
        flux_comparison_output_direct=final_result_obj.flux_comparison_output_direct,
        subset_index_dict=subset_index_dict,
        flux_comparison_xlsx_file_path=final_result_obj.flux_comparison_submitted_output_xlsx_path,
        **other_param_dict)

    metabolic_network_plotting(
        final_solution_data_dict, final_flux_name_index_dict, final_result_obj.metabolic_network_visualization_direct,
        subset_index_dict=subset_index_dict, **other_param_dict)

    output_raw_flux_data(
        final_result_obj.flux_result_output_xlsx_path, final_loss_data_dict,
        final_solution_data_dict, final_flux_name_index_dict,
        final_result_obj.final_information_dict, subset_index_dict=subset_index_dict, other_label_column_dict=None)

    output_predicted_mid_data(
        final_result_obj.mid_prediction_result_output_xlsx_path, final_loss_data_dict, final_predicted_data_dict,
        final_target_experimental_mid_data_dict, final_information_dict, subset_index_dict=subset_index_dict)


def multi_tissue_experimental_mid_prediction(
        experiment_name, complex_predicted_data_dict, final_target_experimental_mid_data_dict,
        solver_dict, mid_name_tissue_name_dict, final_information_dict, result_process_information_dict_analysis,
        tissue_mid_name_list_dict_constructor, prediction_output_direct, subset_index_dict=None,
        predicted_mid_figure_config_dict=None, **kwargs):
    if predicted_mid_figure_config_dict is None:
        predicted_mid_figure_config_dict = {}

    # Extract legend_color_dict if present
    if 'legend_color_dict' in predicted_mid_figure_config_dict:
        legend_color_dict = predicted_mid_figure_config_dict.pop('legend_color_dict')
    else:
        legend_color_dict = None

    # Create and call the template function
    experimental_mid_prediction_plotting_func_template(
        tissue_mid_name_list_dict_constructor,
        result_process_information_dict_analysis,
        final_information_dict, mid_name_tissue_name_dict, solver_dict,
        legend_color_dict=legend_color_dict,
        **predicted_mid_figure_config_dict
    )(
        experiment_name, complex_predicted_data_dict, final_target_experimental_mid_data_dict,
        prediction_output_direct, subset_index_dict=subset_index_dict)


def flux_comparison_parameter_generator(
        experiment_name, solution_array_dict, final_flux_name_index_dict, final_information_dict, flux_name_tissue_name_dict,
        result_process_information_dict_analysis, flux_comparison_output_direct, tissue_flux_name_list_dict_constructor,
        subset_index_dict=None, flux_comparison_figure_config_dict=None, flux_comparison_xlsx_file_path=None,
        **kwargs):
    if flux_comparison_figure_config_dict is None:
        flux_comparison_figure_config_dict = {}

    # Create and call the template function
    organized_mouse_data_dict = flux_comparison_plotting_func_template(
        tissue_flux_name_list_dict_constructor,
        flux_name_tissue_name_dict,
        result_process_information_dict_analysis,
        reversible_flux_pair_dict_generator,
        tissue_specific_reversible_flux_title_constructor_generator,
        **flux_comparison_figure_config_dict
    )(
        experiment_name, solution_array_dict, final_flux_name_index_dict, final_information_dict,
        flux_comparison_output_direct, subset_index_dict=subset_index_dict)

    if flux_comparison_xlsx_file_path is not None:
        output_flux_comparison_data(
            flux_comparison_xlsx_file_path, organized_mouse_data_dict,
            other_label_column_dict=None)


def metabolic_network_plotting(
        final_solution_data_dict, final_flux_name_index_dict, figure_output_direct, subset_index_dict=None,
        metabolic_network_parameter_generator=None, **kwargs):
    return None
    (
        experimental_mid_metabolite_set, experimental_mixed_mid_metabolite_set, biomass_metabolite_set,
        input_metabolite_set, c13_labeling_metabolite_set, boundary_flux_set, infusion
    ) = metabolic_network_parameter_generator()
    for raw_result_label, raw_solution_data_array in final_solution_data_dict.items():
        if subset_index_dict is not None:
            subset_index = subset_index_dict[raw_result_label]
            solution_data_array = raw_solution_data_array[subset_index]
        else:
            solution_data_array = raw_solution_data_array
        current_data_array = solution_data_array.mean(axis=0)
        current_reaction_value_dict = {
            flux_name: current_data_array[flux_index]
            for flux_name, flux_index in final_flux_name_index_dict[raw_result_label].items()}

        output_file_path = f'{figure_output_direct}/metabolic_network_{raw_result_label}.pdf'
        figure_name = f'metabolic_network_{raw_result_label}'
        figure_plotting.metabolic_flux_model_function(
            figure_output_direct, figure_name,
            input_metabolite_set, c13_labeling_metabolite_set, experimental_mid_metabolite_set,
            experimental_mixed_mid_metabolite_set,
            biomass_metabolite_set, boundary_flux_set, current_reaction_value_dict=current_reaction_value_dict,
            infusion=infusion, figure_size=(8.5, 8.5))


def experimental_mid_and_raw_data_plotting(
        experimental_mfa_data_obj_dict, result_information_dict, final_result_obj, result_information_dict_analysis=None,
        tissue_mid_name_list_dict_constructor=None, experimental_figure_config_dict=None):
    if experimental_figure_config_dict is None:
        experimental_figure_config_dict = {}
    experimental_data_output_direct = final_result_obj.experimental_data_output_direct

    experimental_data_plotting_func_template(
        tissue_mid_name_list_dict_constructor, result_information_dict_analysis,
        major_key_file_name_func=lambda major_key: f'target_metabolites_{major_key}_labeling',
        **experimental_figure_config_dict)(
            experimental_mfa_data_obj_dict, result_information_dict, experimental_data_output_direct,
            complete_mid_std_data_vector_dict=None)

    raw_experimental_mid_data_submitted_output_xlsx_path = final_result_obj.raw_experimental_mid_data_submitted_output_xlsx_path
    output_raw_experimental_mid_data(
        raw_experimental_mid_data_submitted_output_xlsx_path, experimental_mfa_data_obj_dict, result_information_dict,
        other_label_row_dict=None)

