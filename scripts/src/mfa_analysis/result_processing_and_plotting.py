from .config import (
    GeneralFinalResult, CoreConstants, multiple_repeat_result_processing, GeneralKeywords,
    mid_name_process, check_if_mix_flux, tissue_name_breakdown, loss_data_distribution_plotting,
    time_distribution_plotting, DataSource, FigureData, FigureDataKeywords,
    tissue_specific_name_constructor, np, reversible_flux_pair_dict_generator, common_flux_comparison_func,
    tissue_specific_reversible_flux_title_constructor_generator, check_and_mkdir_of_direct, figure_plotting,
    Direct, experimental_data_plotting_func_template, TransformDict, output_raw_flux_data, output_predicted_mid_data,
)
from .common_functions import (
    tissue_specific_metabolite_mid_list_constructor, tissue_specific_flux_list_constructor)


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
        # new_target_mid_data_dict = {}
        # for raw_mid_name, mid_data_value in target_mid_data_dict.items():
        #     new_mid_name = raw_mid_name.replace(f'{CoreConstants.mix_ratio_prefix}(', '')
        #     new_mid_name = new_mid_name.replace(f'{CoreConstants.convolution_id}(', '')
        #     new_mid_name = new_mid_name.replace(f')', '')
        #     new_target_mid_data_dict[new_mid_name] = mid_data_value
        # return new_target_mid_data_dict
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

    def load_previous_results(self, raw_result_label):
        return self._load_previous_results(raw_result_label, out_of_order=True)


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
        subset_index_dict=subset_index_dict, **other_param_dict)

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
    tissue_mid_name_list_dict, display_mid_name_dict, display_tissue_name_dict = tissue_mid_name_list_dict_constructor(
        result_process=True)

    formatted_mean_data_dict = {}
    formatted_std_data_dict = {}
    data_len_dict = {}
    subplot_name_dict = {}
    color_dict = {}
    default_experimental_transparency = 0.3
    result_information_dict_list = result_process_information_dict_analysis(
        final_information_dict, target=GeneralKeywords.mid_prediction)
    # assert len(result_information_dict_list) == len(complex_predicted_data_dict)
    # Filter out information for results that are missing in the predicted data (e.g. failed runs with 0 solutions)
    result_information_dict_list = [
        item for item in result_information_dict_list
        if item[0] in complex_predicted_data_dict
    ]

    for (
            result_label, major_key, minor_key_list, minor_key_str, current_color, order_index
    ) in result_information_dict_list:
        result_specific_predicted_data_dict = complex_predicted_data_dict[result_label]
        current_raw_emu_name_experimental_name_dict = solver_dict[result_label].emu_name_experimental_name_dict
        current_emu_name_experimental_name_dict = {
            mid_name_process(raw_emu_name): experiment_name
            for raw_emu_name, experiment_name in current_raw_emu_name_experimental_name_dict.items()}
        current_solver_target_experimental_mid_dict = {
            mid_name_process(raw_mid_name): mid_vector
            for raw_mid_name, mid_vector in solver_dict[result_label].target_experimental_mid_data_dict.items()
        }
        predicted_key_str = minor_key_str
        experimental_key_str = f'{minor_key_str}_{GeneralKeywords.experimental}'
        color_dict.update({
            predicted_key_str: current_color,
            experimental_key_str: current_color.transparency_mix(default_experimental_transparency),
        })
        current_target_experimental_mid_data_dict = final_target_experimental_mid_data_dict[result_label]
        # result_label_optimization_state_dict.update({
        #     result_label: (result_label, Keywords.optimized),
        #     experimental_result_label: (result_label, Keywords.experimental),
        # })
        # TODO: Need to mark the existing but manually excluded MFA data, similar with experimental_data_plotting.
        for mid_title, current_predicted_data_array_list in result_specific_predicted_data_dict.items():
            if mid_title in mid_name_tissue_name_dict:
                tissue_name, raw_mid_name = mid_name_tissue_name_dict[mid_title]
            else:
                tissue_name, raw_mid_name = tissue_name_breakdown(mid_title)
                mid_name_tissue_name_dict[mid_title] = (tissue_name, raw_mid_name)
            current_experimental_name = current_emu_name_experimental_name_dict[mid_title]
            updated_tissue_specific_mid_name = tissue_specific_name_constructor(current_experimental_name, tissue_name)
            current_experimental_mid = current_solver_target_experimental_mid_dict[mid_title]
            if major_key not in formatted_mean_data_dict:
                formatted_mean_data_dict[major_key] = {}
                formatted_std_data_dict[major_key] = {}
            current_formatted_mean_data_dict = formatted_mean_data_dict[major_key]
            current_formatted_std_data_dict = formatted_std_data_dict[major_key]
            if updated_tissue_specific_mid_name not in current_formatted_mean_data_dict:
                current_formatted_mean_data_dict[updated_tissue_specific_mid_name] = {}
                current_formatted_std_data_dict[updated_tissue_specific_mid_name] = {}
                if display_mid_name_dict is not None and raw_mid_name in display_mid_name_dict:
                    display_mid_name = display_mid_name_dict[raw_mid_name]
                else:
                    display_mid_name = current_experimental_name
                display_tissue_mid_name = f'{tissue_name}: {display_mid_name}'
                subplot_name_dict[updated_tissue_specific_mid_name] = display_tissue_mid_name
            current_predicted_data_array = np.array(current_predicted_data_array_list)
            if subset_index_dict is not None:
                target_predicted_data_array = current_predicted_data_array[subset_index_dict[result_label], :]
            else:
                target_predicted_data_array = current_predicted_data_array
            if current_experimental_name not in data_len_dict:
                data_len_dict[current_experimental_name] = current_predicted_data_array.shape[1]
            current_formatted_mean_data_dict[updated_tissue_specific_mid_name][predicted_key_str] = (
                target_predicted_data_array.mean(axis=0))
            current_formatted_std_data_dict[updated_tissue_specific_mid_name][predicted_key_str] = (
                target_predicted_data_array.std(axis=0))
            current_formatted_mean_data_dict[updated_tissue_specific_mid_name][experimental_key_str] = (
                current_experimental_mid
                # final_target_experimental_mid_data_dict[result_label][mid_title]
            )

    mean_data_dict = {}
    std_data_dict = {}
    mid_name_list_dict = {}
    if 'legend_color_dict' in predicted_mid_figure_config_dict:
        legend_color_dict = predicted_mid_figure_config_dict.pop('legend_color_dict')
    else:
        legend_color_dict = None
    for major_key, current_formatted_mean_data in formatted_mean_data_dict.items():
        current_formatted_std_data = formatted_std_data_dict[major_key]
        current_mid_name_list, current_mean_data_dict, current_std_data_dict = tissue_specific_metabolite_mid_list_constructor(
            current_formatted_mean_data, current_formatted_std_data, data_len_dict, tissue_mid_name_list_dict)
        mean_data_dict[major_key] = current_formatted_mean_data
        std_data_dict[major_key] = current_formatted_std_data
        mid_name_list_dict[major_key] = current_mid_name_list

        target_row_num = len(current_mid_name_list)
        target_col_num = max([len(each_row_mid_name_list) for each_row_mid_name_list in current_mid_name_list])
        current_legend_color_dict = legend_color_dict[major_key] if legend_color_dict is not None else None

        figure_plotting.mid_comparison_plotting(
            figure_output_direct=prediction_output_direct, figure_name=f'predicted_mid_comparison_{major_key}',
            complete_data_dict=current_mean_data_dict, target_emu_name_nested_list=current_mid_name_list,
            target_row_num=target_row_num, target_col_num=target_col_num, error_bar_data_dict=current_std_data_dict,
            color_dict=color_dict, legend_color_dict=current_legend_color_dict, y_lim=(0, 1),
            **predicted_mid_figure_config_dict)

    figure_raw_data = FigureData(FigureDataKeywords.mid_comparison, experiment_name)
    figure_raw_data.save_data(
        mean_data_dict=mean_data_dict,
        std_data_dict=std_data_dict,
        mid_name_list=mid_name_list_dict,
        subplot_name_dict=subplot_name_dict,
    )


def flux_comparison_parameter_generator(
        experiment_name, solution_array_dict, final_flux_name_index_dict, final_information_dict, flux_name_tissue_name_dict,
        result_process_information_dict_analysis, flux_comparison_output_direct, tissue_flux_name_list_dict_constructor,
        subset_index_dict=None, flux_comparison_figure_config_dict=None, **kwargs):
    if flux_comparison_figure_config_dict is None:
        flux_comparison_figure_config_dict = {}

    (
        all_tissue_flux_name_list_dict, display_flux_dict, specific_tissue_reversible_flux_pair_dict
    ) = tissue_flux_name_list_dict_constructor()

    preprocessed_flux_data_dict = {}
    for result_label, current_solution_array in solution_array_dict.items():
        current_flux_name_index_dict = final_flux_name_index_dict[result_label]
        reversible_flux_pair_dict, all_flux_name_list = reversible_flux_pair_dict_generator(
            current_flux_name_index_dict, flux_name_tissue_name_dict, specific_tissue_reversible_flux_pair_dict)
        if subset_index_dict is not None:
            target_solution_array = current_solution_array[subset_index_dict[result_label], :]
        else:
            target_solution_array = current_solution_array
        preprocessed_flux_data_dict[result_label] = (
            current_flux_name_index_dict, all_flux_name_list, target_solution_array)

    organized_mouse_data_dict = {}
    color_dict = {}
    name_dict = {}
    for (
            result_label, major_key, minor_key_list, minor_key_str, current_color, order_index
    ) in result_process_information_dict_analysis(final_information_dict, target=GeneralKeywords.flux_comparison):
        (
            current_flux_name_index_dict, all_flux_name_list, target_solution_array
        ) = preprocessed_flux_data_dict[result_label]
        if major_key not in organized_mouse_data_dict:
            organized_mouse_data_dict[major_key] = {}
            color_dict[major_key] = {}
            name_dict[major_key] = {}
        color_dict[major_key][minor_key_str] = current_color
        name_dict[major_key][minor_key_str] = minor_key_str
        common_flux_comparison_func(
            all_flux_name_list, current_flux_name_index_dict, target_solution_array,
            organized_mouse_data_dict[major_key], minor_key_str, extra_group_key_name='',
            reversible_flux_title_constructor=tissue_specific_reversible_flux_title_constructor_generator(
                flux_name_tissue_name_dict))

    # for result_label, current_solution_array in solution_array_dict.items():
    #     current_flux_name_index_dict = flux_name_index_dict[result_label]
    #     reversible_flux_pair_dict, all_flux_name_list = reversible_flux_pair_dict_generator(
    #         current_flux_name_index_dict, flux_name_tissue_name_dict, specific_tissue_reversible_flux_pair_dict)
    #     major_key, minor_key_list, minor_key_str, current_color, order_index = result_process_information_dict_analysis(
    #         information_dict[result_label])
    #     if major_key not in organized_mouse_data_dict:
    #         organized_mouse_data_dict[major_key] = {}
    #         color_dict[major_key] = {}
    #         name_dict[major_key] = {}
    #     color_dict[major_key][minor_key_str] = current_color
    #     name_dict[major_key][minor_key_str] = minor_key_str
    #     if subset_index_dict is not None:
    #         target_solution_array = current_solution_array[subset_index_dict[result_label], :]
    #     else:
    #         target_solution_array = current_solution_array
    #     common_flux_comparison_func(
    #         all_flux_name_list, current_flux_name_index_dict, target_solution_array,
    #         organized_mouse_data_dict[major_key], minor_key_str, extra_group_key_name='',
    #         reversible_flux_title_constructor=tissue_specific_reversible_flux_title_constructor_generator(
    #             flux_name_tissue_name_dict))

    for major_key, current_data_dict in organized_mouse_data_dict.items():
        current_flux_name_nested_list, tissue_specific_display_flux_dict = tissue_specific_flux_list_constructor(
            all_tissue_flux_name_list_dict, display_flux_dict)
        figure_plotting.distribution_comparison_plotting(
            figure_output_direct=flux_comparison_output_direct, figure_name=f'flux_comparison_{major_key}',
            complete_data_dict=current_data_dict, flux_name_nested_list=current_flux_name_nested_list,
            display_flux_dict=tissue_specific_display_flux_dict, name_dict=name_dict[major_key],
            color_dict=color_dict[major_key], **flux_comparison_figure_config_dict)

    figure_raw_data = FigureData(FigureDataKeywords.flux_comparison, experiment_name)
    figure_raw_data.save_data(
        organized_diet_mouse_data_dict=organized_mouse_data_dict,)


def mfa_data_dict_preprocess(
        mouse_id_mfa_data_dict, std_mid_data_vector_dict, display_mid_name_dict, display_tissue_name_dict,
        specific_dimension_metabolite_dict=None):
    if specific_dimension_metabolite_dict is None:
        specific_dimension_metabolite_dict = {}
    formatted_mean_data_dict = {}
    formatted_std_data_dict = {}
    subplot_name_dict = {}
    data_len_dict = {
        metabolite_name: dimension + 1
        for metabolite_name, dimension in specific_dimension_metabolite_dict.items()}
    for mouse_id, mfa_data_obj in mouse_id_mfa_data_dict.items():
        for tissue_metabolite_name, mid_data_obj in mfa_data_obj.experimental_mid_data_obj_dict.items():
            metabolite_name = mid_data_obj.name
            mid_data_vector = mid_data_obj.data_vector
            tissue = mid_data_obj.tissue
            if isinstance(tissue, (tuple, list)) and len(tissue) == 1:
                tissue = tissue[0]
            if metabolite_name not in data_len_dict:
                data_len_dict[metabolite_name] = len(mid_data_vector)
            new_tissue_specific_mid_name = tissue_specific_name_constructor(metabolite_name, tissue)
            if new_tissue_specific_mid_name not in formatted_mean_data_dict:
                formatted_mean_data_dict[new_tissue_specific_mid_name] = {}
                formatted_std_data_dict[new_tissue_specific_mid_name] = {}
            formatted_mean_data_dict[new_tissue_specific_mid_name][mouse_id] = mid_data_vector
            try:
                current_std_data_vector = std_mid_data_vector_dict[mouse_id][tissue_metabolite_name]
            except KeyError:
                current_std_data_vector = None
            formatted_std_data_dict[new_tissue_specific_mid_name][mouse_id] = current_std_data_vector
            if new_tissue_specific_mid_name not in subplot_name_dict:
                display_mid_name = (
                    display_mid_name_dict[metabolite_name] if metabolite_name in display_mid_name_dict
                    else metabolite_name)
                if tissue is None:
                    display_tissue_mid_name = display_mid_name
                else:
                    display_tissue_name = (
                        display_tissue_name_dict[tissue] if tissue in display_tissue_name_dict
                        else tissue)
                    display_tissue_mid_name = f'{display_tissue_name}: {display_mid_name}'
                subplot_name_dict[new_tissue_specific_mid_name] = display_tissue_mid_name
    return formatted_mean_data_dict, formatted_std_data_dict, subplot_name_dict, data_len_dict


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
        mouse_id_mfa_data_dict, result_information_dict, final_result_obj, result_information_dict_analysis=None,
        tissue_mid_name_list_dict_constructor=None, experimental_figure_config_dict=None, average_data=False):
    if experimental_figure_config_dict is None:
        experimental_figure_config_dict = {}
    experimental_data_output_direct = final_result_obj.experimental_data_output_direct

    experimental_data_plotting_func_template(
        mouse_id_mfa_data_dict, tissue_mid_name_list_dict_constructor, result_information_dict_analysis,
        major_key_file_name_func=lambda major_key: f'target_metabolites_{major_key}_labeling',
        **experimental_figure_config_dict)(
            mouse_id_mfa_data_dict, result_information_dict, experimental_data_output_direct,
            complete_mid_std_data_vector_dict=None)


def plotting_mid_data(
        data_dict, error_bar_data_dict=None, mouse_id_diet_dict=None, mid_name_list=None, subplot_name_dict=None,
        low_height=False, figure_title=None, figure_name=None, figure_output_direct=None,
        color_dict=None, name_dict=None, diet_data_dict=False):
    figure_size = (8.5, 11)
    ParameterName = figure_plotting.ParameterName

    if mid_name_list is None:
        raise ValueError()
    diet_mouse_id_list_dict = {}
    for mouse_id, diet_name in mouse_id_diet_dict.items():
        if diet_name not in diet_mouse_id_list_dict:
            diet_mouse_id_list_dict[diet_name] = []
        diet_mouse_id_list_dict[diet_name].append(mouse_id)
    if diet_data_dict:
        final_color_dict = color_dict
    else:
        mouse_id_color_dict = {}
        for diet_name, diet_color in color_dict.items():
            for mouse_id in diet_mouse_id_list_dict[diet_name]:
                mouse_id_color_dict[mouse_id] = diet_color
        final_color_dict = mouse_id_color_dict
    # name_dict = {mouse_id: mouse_id for mouse_id in mouse_id_diet_dict.keys()}
    legend_color_dict = color_dict
    legend_name_dict = name_dict
    if subplot_name_dict is None:
        subplot_name_dict = {}
    figure_data_parameter_dict = {
        ParameterName.figure_data: (data_dict, error_bar_data_dict),
        ParameterName.figure_title: figure_title,
        ParameterName.name_dict: legend_name_dict,
        ParameterName.legend_color_dict: legend_color_dict,
        ParameterName.color_dict: final_color_dict,
        ParameterName.mid_name_list: mid_name_list,
        ParameterName.subplot_name_dict: subplot_name_dict,
        ParameterName.low_height: low_height,
        ParameterName.figure_config_dict: {},
    }
    figure_plotting.mid_prediction_function(
        figure_output_direct=figure_output_direct, figure_name=figure_name,
        figure_data_parameter_dict=figure_data_parameter_dict, figure_size=figure_size)


