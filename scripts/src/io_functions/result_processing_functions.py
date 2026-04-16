from .config import (
    check_and_mkdir_of_direct, np, Keywords, FigureData, FigureDataKeywords, mid_name_process,
    tissue_specific_name_constructor, TransformDict, CoreConstants)
from .plotting_functions import figure_plotting
from .common_functions import (
    tissue_specific_metabolite_mid_list_constructor)


def add_empty_obj(root_dict, final_empty_class, *args):
    current_dict = root_dict
    total_arg_num = len(args)
    for label_index, label in enumerate(args):
        if label not in current_dict:
            if label_index == total_arg_num - 1:
                current_dict[label] = final_empty_class()
            else:
                current_dict[label] = {}
        current_dict = current_dict[label]
    return current_dict


def sort_dict_and_remove_sorting_index(raw_dict):
    ordered_list = list(raw_dict.items())
    ordered_list.sort(key=lambda x: x[1][0])
    new_dict = {}
    for key, (value_index, value) in ordered_list:
        new_dict[key] = value
    return new_dict


def convert_tmp_data_to_final_data_dict(tmp_data_dict):
    final_data_dict = {}
    for major_key, each_major_key_tmp_data_dict in tmp_data_dict.items():
        final_data_dict[major_key] = {}
        for metabolite_name, each_metabolite_data_dict in each_major_key_tmp_data_dict.items():
            final_data_dict[major_key][metabolite_name] = sort_dict_and_remove_sorting_index(
                each_metabolite_data_dict)
    return final_data_dict


def experimental_data_plotting_func_template(
        tissue_mid_name_list_dict_constructor, major_minor_key_analysis_func,
        major_key_file_name_func=lambda major_key: f'target_metabolites_{major_key}_labeling',
        legend_color_dict=None, **figure_param_dict):
    def preprocess_all_result_label_data(
            complete_experimental_mfa_data_obj_dict, complete_mid_std_data_vector_dict):
        mid_data_dict = {}
        raw_data_dict = {}
        std_data_dict = {}
        mfa_exclusion_dict = {}
        for result_label, experimental_mfa_data_obj in complete_experimental_mfa_data_obj_dict.items():
            if result_label not in mid_data_dict:
                mid_data_dict[result_label] = {}
                raw_data_dict[result_label] = {}
                std_data_dict[result_label] = {}
                mfa_exclusion_dict[result_label] = {}
            for metabolite_name, mid_data_obj in experimental_mfa_data_obj.experimental_mid_data_obj_dict.items():
                mid_data_dict[result_label][metabolite_name] = mid_data_obj.data_vector
                raw_data_dict[result_label][metabolite_name] = mid_data_obj.raw_data_vector
                if complete_mid_std_data_vector_dict is not None:
                    try:
                        current_std_vector = complete_mid_std_data_vector_dict[result_label][metabolite_name]
                    except KeyError:
                        current_std_vector = None
                    std_data_dict[result_label][metabolite_name] = current_std_vector
                mfa_exclusion_dict[result_label][metabolite_name] = mid_data_obj.excluded_from_mfa
        return mid_data_dict, raw_data_dict, std_data_dict, mfa_exclusion_dict

    def experimental_data_plotting(
            complete_experimental_mfa_data_obj_dict, complete_result_information_dict, output_direct,
            complete_mid_std_data_vector_dict=None,):
        (
            preprocess_mid_data_dict, preprocess_raw_data_dict, preprocess_std_data_dict, preprocess_mfa_exclusion_dict
        ) = preprocess_all_result_label_data(
            complete_experimental_mfa_data_obj_dict, complete_mid_std_data_vector_dict)
        tmp_mid_data_dict = {}
        tmp_raw_data_dict = {}
        tmp_std_data_dict = {}
        tmp_color_dict = {}
        tmp_mfa_excluded_dict = {}
        for (
                result_label, major_key, minor_key_list, minor_key_str, current_color, order_index
        ) in major_minor_key_analysis_func(complete_result_information_dict):
            for metabolite_name, current_metabolite_mid_data in preprocess_mid_data_dict[result_label].items():
                current_tmp_mid_data_dict = add_empty_obj(tmp_mid_data_dict, dict, major_key, metabolite_name)
                current_tmp_raw_data_dict = add_empty_obj(tmp_raw_data_dict, dict, major_key, metabolite_name)
                current_tmp_mid_data_dict[minor_key_str] = (order_index, current_metabolite_mid_data)
                current_tmp_raw_data_dict[minor_key_str] = (order_index, preprocess_raw_data_dict[result_label][metabolite_name])
                if result_label in tmp_std_data_dict:
                    current_tmp_std_data_dict = add_empty_obj(tmp_std_data_dict, dict, major_key, metabolite_name)
                    current_tmp_std_data_dict[minor_key_str] = (order_index, preprocess_std_data_dict[result_label][metabolite_name])
                current_mfa_excluded_dict = add_empty_obj(tmp_mfa_excluded_dict, dict, major_key, metabolite_name)
                current_mfa_excluded_dict[minor_key_str] = (order_index, preprocess_mfa_exclusion_dict[result_label][metabolite_name])
                current_tmp_color_dict = add_empty_obj(tmp_color_dict, dict, major_key)
                if minor_key_str not in current_tmp_color_dict:
                    current_tmp_color_dict[minor_key_str] = (order_index, current_color)
        mid_data_dict_for_plotting = convert_tmp_data_to_final_data_dict(tmp_mid_data_dict)
        raw_data_dict_for_plotting = convert_tmp_data_to_final_data_dict(tmp_raw_data_dict)
        std_data_dict_for_plotting = convert_tmp_data_to_final_data_dict(tmp_std_data_dict)
        mfa_excluded_dict_for_plotting = convert_tmp_data_to_final_data_dict(tmp_mfa_excluded_dict)
        complete_color_dict = {}
        for major_key, each_major_key_tmp_color_dict in tmp_color_dict.items():
            complete_color_dict[major_key] = sort_dict_and_remove_sorting_index(each_major_key_tmp_color_dict)

        (
            tissue_mid_name_list_dict, specific_dimension_metabolite_dict,
            display_mid_name_dict, display_tissue_name_dict) = tissue_mid_name_list_dict_constructor()
        # target_row_num = len(target_emu_name_nested_list)
        # target_col_num = len(target_emu_name_nested_list[0])
        for major_key, each_labeling_mid_data_dict_for_plotting in mid_data_dict_for_plotting.items():
            each_labeling_mfa_excluded_dict_for_plotting = mfa_excluded_dict_for_plotting[major_key]
            (
                target_emu_name_nested_list, title_dict, _, mfa_excluded_text_list, mfa_excluded_text_loc_list
            ) = mfa_data_dict_preprocess(
                each_labeling_mid_data_dict_for_plotting, display_mid_name_dict, display_tissue_name_dict, tissue_mid_name_list_dict,
                each_labeling_mfa_excluded_dict_for_plotting,
                specific_dimension_metabolite_dict=specific_dimension_metabolite_dict, **figure_param_dict)
            target_row_num = len(target_emu_name_nested_list)
            target_col_num = len(target_emu_name_nested_list[0])

            each_labeling_raw_data_dict_for_plotting = raw_data_dict_for_plotting[major_key]
            color_dict = complete_color_dict[major_key]
            if legend_color_dict is not None:
                current_legend_color_dict = legend_color_dict[major_key]
            else:
                current_legend_color_dict = None
            for raw_data in (False, True):
                if raw_data:
                    parent_direct = 'raw_data'
                    complete_data_dict = each_labeling_raw_data_dict_for_plotting
                    ylim = (0, None)
                    supplementary_text_list = None
                    supplementary_text_loc_list = None
                else:
                    parent_direct = 'mid_data'
                    complete_data_dict = each_labeling_mid_data_dict_for_plotting
                    ylim = (0, 1)
                    supplementary_text_list = mfa_excluded_text_list
                    supplementary_text_loc_list = mfa_excluded_text_loc_list
                current_title = major_key_file_name_func(major_key)
                current_cell_line_output_direct = '{}/{}'.format(output_direct, parent_direct)
                check_and_mkdir_of_direct(current_cell_line_output_direct)
                if len(std_data_dict_for_plotting) == 0:
                    std_data_dict = None
                else:
                    std_data_dict = std_data_dict_for_plotting
                figure_plotting.mid_comparison_plotting(
                    figure_output_direct=current_cell_line_output_direct, figure_name=current_title,
                    complete_data_dict=complete_data_dict, target_emu_name_nested_list=target_emu_name_nested_list,
                    target_row_num=target_row_num, target_col_num=target_col_num,
                    error_bar_data_dict=std_data_dict, color_dict=color_dict, y_lim=ylim,
                    legend_color_dict=current_legend_color_dict, title_dict=title_dict,
                    supplementary_text_list=supplementary_text_list,
                    supplementary_text_loc_list=supplementary_text_loc_list, **figure_param_dict)

    return experimental_data_plotting


def mfa_data_dict_preprocess(
        formatted_raw_data_dict, display_mid_name_dict, display_tissue_name_dict, tissue_mid_name_list_dict,
        mfa_excluded_dict, specific_dimension_metabolite_dict=None, mfa_excluded_text_loc=(0.5, 0.5), **kwargs):
    if specific_dimension_metabolite_dict is None:
        specific_dimension_metabolite_dict = {}
    subplot_name_dict = {}
    data_len_dict = {
        metabolite_name: dimension + 1
        for metabolite_name, dimension in specific_dimension_metabolite_dict.items()}
    specific_tissue_sep = CoreConstants.specific_tissue_sep
    for new_tissue_specific_mid_name in formatted_raw_data_dict.keys():
        if specific_tissue_sep in new_tissue_specific_mid_name:
            tissue, metabolite_name = new_tissue_specific_mid_name.split(specific_tissue_sep)
        else:
            tissue = None
            metabolite_name = new_tissue_specific_mid_name
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
    mid_name_list = []
    mfa_excluded_text_list = []
    mfa_excluded_text_loc_list = []
    filtered_raw_data_dict = {}
    for tissue_name, each_tissue_mid_name_list in tissue_mid_name_list_dict.items():
        for mid_name_row in each_tissue_mid_name_list:
            new_mid_name_row = []
            for metabolite_name in mid_name_row:
                new_tissue_specific_mid_name = tissue_specific_name_constructor(metabolite_name, tissue_name)
                try:
                    this_mid_data_dict = formatted_raw_data_dict[new_tissue_specific_mid_name]
                except KeyError:
                    this_mid_data_dict = None
                    try:
                        mid_data_len = data_len_dict[metabolite_name]
                    except KeyError:
                        mid_data_len = None
                    new_mid_name_row.append(mid_data_len)
                    current_mfa_excluded_text = None
                else:
                    new_mid_name_row.append(new_tissue_specific_mid_name)
                    if new_tissue_specific_mid_name in mfa_excluded_dict:
                        this_mid_excluded_dict = mfa_excluded_dict[new_tissue_specific_mid_name]
                        rm_minor_key_list = [
                            minor_key for minor_key, excluded in this_mid_excluded_dict.items() if excluded]
                    else:
                        rm_minor_key_list = []
                    if len(rm_minor_key_list) > 0:
                        current_mfa_excluded_text = f'Excluded from MFA:\n{"\n".join(rm_minor_key_list)}'
                    else:
                        current_mfa_excluded_text = None
                mfa_excluded_text_list.append(current_mfa_excluded_text)
                mfa_excluded_text_loc_list.append(mfa_excluded_text_loc)
                filtered_raw_data_dict[new_tissue_specific_mid_name] = this_mid_data_dict
            mid_name_list.append(new_mid_name_row)
    title_dict = TransformDict(subplot_name_dict)
    return (
        mid_name_list, title_dict, filtered_raw_data_dict, mfa_excluded_text_list, mfa_excluded_text_loc_list)


def mid_name_list_generator_for_multiple_labeling_substrate(raw_metabolite_list, labeling_substrate_name_list):
    final_metabolite_list = []
    for labeling_substrate_name in labeling_substrate_name_list:
        for raw_metabolite_row in raw_metabolite_list:
            combined_metabolite_row = []
            for raw_metabolite_name in raw_metabolite_row:
                if raw_metabolite_name is None:
                    combined_metabolite_name = None
                else:
                    combined_metabolite_name = f'{labeling_substrate_name}_{raw_metabolite_name}'
                combined_metabolite_row.append(combined_metabolite_name)
            final_metabolite_list.append(combined_metabolite_row)
    return final_metabolite_list


def loss_data_distribution_plotting(
        experiment_name, loss_data_dict, output_direct=None, loss_percentile=None, select_num=None):
    if output_direct is not None:
        # group_violin_box_distribution_plot(
        #     {'loss_distribution': loss_data_dict}, nested_color_dict=None,
        #     nested_median_color_dict=None, cutoff_dict=None, title_dict=None,
        #     output_direct=output_direct, ylim=None, xaxis_rotate=True,
        #     figsize=None, figure_type='box')
        figure_plotting.violin_box_distribution_function(
            figure_output_direct=output_direct, figure_name='loss_distribution',
            complete_data_dict={'loss_distribution': loss_data_dict}, figure_type='box')
    if loss_percentile is not None or select_num is not None:
        if loss_percentile is not None and select_num is not None:
            raise ValueError()
        subset_index_dict = {}
        filtered_loss_data_dict = {}
        for result_label, loss_array in loss_data_dict.items():
            filtered_index = select_solutions(loss_array, loss_percentile, select_num)
            subset_index_dict[result_label] = filtered_index
            filtered_loss_data_dict[result_label] = loss_array[filtered_index]
        if output_direct is not None:
            # group_violin_box_distribution_plot(
            #     {'filtered_loss_distribution': filtered_loss_data_dict}, nested_color_dict=None,
            #     nested_median_color_dict=None, cutoff_dict=None, title_dict=None,
            #     output_direct=output_direct, ylim=None, xaxis_rotate=True,
            #     figsize=None, figure_type='box')
            figure_plotting.violin_box_distribution_function(
                figure_output_direct=output_direct, figure_name='filtered_loss_distribution',
                complete_data_dict={'filtered_loss_distribution': filtered_loss_data_dict}, figure_type='box')
    else:
        subset_index_dict = None
        filtered_loss_data_dict = loss_data_dict
    if output_direct is not None:
        figure_raw_data = FigureData(FigureDataKeywords.loss_data_comparison, experiment_name)
        figure_raw_data.save_data(
            loss_data_dict=loss_data_dict,
            filtered_loss_data_dict=filtered_loss_data_dict)
    return subset_index_dict


def time_distribution_plotting(
        experiment_name, time_data_dict, output_direct=None, subset_index_dict=None, parallel_num=None):
    filtered_time_data_dict = {}
    if parallel_num is not None:
        for result_label, time_data_array in time_data_dict.items():
            time_data_dict[result_label] = time_data_array / parallel_num
    if output_direct is not None:
        # group_violin_box_distribution_plot(
        #     {'time_distribution': time_data_dict}, nested_color_dict=None,
        #     nested_median_color_dict=None, cutoff_dict=None, title_dict=None,
        #     output_direct=output_direct, ylim=None, xaxis_rotate=True,
        #     figsize=None, figure_type='box')
        figure_plotting.violin_box_distribution_function(
            figure_output_direct=output_direct, figure_name='time_distribution',
            complete_data_dict={'time_distribution': time_data_dict}, figure_type='box')
        if subset_index_dict is not None:
            for result_label, time_array in time_data_dict.items():
                filtered_index = subset_index_dict[result_label]
                filtered_time_data_dict[result_label] = time_array[filtered_index]
            # group_violin_box_distribution_plot(
            #     {'filtered_time_distribution': filtered_time_data_dict}, nested_color_dict=None,
            #     nested_median_color_dict=None, cutoff_dict=None, title_dict=None,
            #     output_direct=output_direct, ylim=None, xaxis_rotate=True,
            #     figsize=None, figure_type='box')
            figure_plotting.violin_box_distribution_function(
                figure_output_direct=output_direct, figure_name='filtered_time_distribution',
                complete_data_dict={'filtered_time_distribution': filtered_time_data_dict}, figure_type='box')
    if output_direct is not None:
        figure_raw_data = FigureData(FigureDataKeywords.time_data_distribution, experiment_name)
        figure_raw_data.save_data(
            time_data_dict=time_data_dict, filtered_time_data_dict=filtered_time_data_dict)


def reconstruct_and_filter_data_dict(
        final_solution_data_dict, final_flux_name_index_dict, final_mapping_dict, subset_index_dict=None):
    def decouple_result_label_tuple(label_tuple, data_dict, data_array):
        current_label = label_tuple[0]
        if len(label_tuple) == 1:
            data_dict[current_label] = data_array
        else:
            if current_label not in data_dict:
                data_dict[current_label] = {}
            decouple_result_label_tuple(label_tuple[1:], data_dict[current_label], data_array)

    reconstructed_solution_data_dict = {}
    reconstructed_flux_name_index_dict = {}
    common_flux_name_index_dict = None
    for raw_result_label, raw_solution_data_array in final_solution_data_dict.items():
        if subset_index_dict is not None:
            subset_index = subset_index_dict[raw_result_label]
            solution_data_array = raw_solution_data_array[subset_index]
        else:
            solution_data_array = raw_solution_data_array
        # if common_flux_name_index_dict is None:
        #     common_flux_name_index_dict = final_flux_name_index_dict[raw_result_label]
        current_flux_name_index_dict = final_flux_name_index_dict[raw_result_label]
        complete_result_label_tuple = final_mapping_dict[raw_result_label]
        decouple_result_label_tuple(complete_result_label_tuple, reconstructed_solution_data_dict, solution_data_array)
        decouple_result_label_tuple(
            complete_result_label_tuple, reconstructed_flux_name_index_dict, current_flux_name_index_dict)
    return reconstructed_solution_data_dict, reconstructed_flux_name_index_dict


def select_solutions(loss_array, loss_percentile=None, select_num=None, index_start=0):
    index_array = np.argsort(loss_array)
    total_num = loss_array.shape[0]
    target_num = total_num
    if loss_percentile is not None:
        target_num = int(loss_percentile * total_num + 0.9999)
    elif select_num is not None:
        target_num = select_num
    filtered_index = index_array[:target_num]
    return filtered_index + index_start


def multiple_repeat_result_processing(
        solver_dict, final_solution_data_dict, final_loss_data_dict, final_predicted_mid_data_dict,
        repeat_division_num=1, each_optimization_num=None, loss_percentile=None, select_num=None):
    def index_start_end_generator(total_num, batch_num, optimization_num, start_index=0):
        if total_num == 0:
            return
        if optimization_num is not None:
            base_num = optimization_num
            if batch_num is None:
                batch_num = total_num // base_num
        else:
            if batch_num is not None and batch_num > total_num:
                batch_num = total_num
            base_num = total_num // batch_num

        if base_num == 0:
            return


        remainder = total_num % base_num
        size_list = [base_num] * batch_num
        current_start_index = start_index
        for current_size in size_list:
            current_end_index = current_start_index + current_size
            yield current_start_index, current_end_index
            current_start_index = current_end_index

    new_loss_data_dict = {}
    new_final_solution_data_dict = {}
    new_final_predicted_mid_data_dict = {}

    for result_label, loss_array in final_loss_data_dict.items():
        if isinstance(repeat_division_num, dict):
            this_result_repeat_division = repeat_division_num[result_label]
        elif isinstance(repeat_division_num, int):
            this_result_repeat_division = repeat_division_num
        elif repeat_division_num is None:
            this_result_repeat_division = None
        else:
            raise TypeError
        if isinstance(each_optimization_num, dict):
            this_result_optimization_num = each_optimization_num[result_label]
        elif isinstance(each_optimization_num, int):
            this_result_optimization_num = each_optimization_num
        elif each_optimization_num is None:
            this_result_optimization_num = None
        else:
            raise TypeError
        solver_obj = solver_dict[result_label]
        solution_array = final_solution_data_dict[result_label]
        predicted_mid_data_dict = final_predicted_mid_data_dict[result_label]
        new_loss_list = []
        new_solution_list = []
        new_mid_data_dict = {}
        total_solution_num = len(loss_array)
        for index_start, index_end in index_start_end_generator(
                total_solution_num, this_result_repeat_division, this_result_optimization_num):
            this_repeat_loss_array = loss_array[index_start:index_end]
            filtered_index = select_solutions(this_repeat_loss_array, loss_percentile, select_num, index_start)
            mean_solution_array = np.mean(solution_array[filtered_index, :], 0)
            calculated_loss = solver_obj.obj_eval(mean_solution_array)
            calculated_mid_data_dict = solver_obj.predict(mean_solution_array)
            new_solution_list.append(mean_solution_array)
            new_loss_list.append(calculated_loss)
            for mid_title, predicted_mid_array in calculated_mid_data_dict.items():
                processed_mid_title = mid_name_process(mid_title)
                if processed_mid_title not in new_mid_data_dict:
                    new_mid_data_dict[processed_mid_title] = []
                new_mid_data_dict[processed_mid_title].append(predicted_mid_array)
        new_loss_data_dict[result_label] = np.array(new_loss_list)
        new_final_solution_data_dict[result_label] = np.array(new_solution_list)
        new_final_predicted_mid_data_dict[result_label] = new_mid_data_dict

    return new_final_solution_data_dict, new_loss_data_dict, new_final_predicted_mid_data_dict


def multiple_repeat_result_processing_new(
        solver_dict, final_solution_data_dict, final_loss_data_dict, final_predicted_mid_data_dict,
        repeat_division_num=1, each_optimization_num=None, loss_percentile=None, select_num=None):
    def index_start_end_generator(total_num, batch_num, optimization_num, start_index=0):
        if total_num == 0:
            return
        if optimization_num is not None:
            base_num = optimization_num
            if batch_num is None:
                batch_num = total_num // base_num
        else:
            if batch_num is not None and batch_num > total_num:
                batch_num = total_num
            base_num = total_num // batch_num

        if base_num == 0:
            return


        remainder = total_num % base_num
        size_list = [base_num] * batch_num
        current_start_index = start_index
        for current_size in size_list:
            current_end_index = current_start_index + current_size
            yield current_start_index, current_end_index
            current_start_index = current_end_index

    new_loss_data_dict = {}
    new_final_solution_data_dict = {}
    new_final_predicted_mid_data_dict = {}

    for result_label, loss_array in final_loss_data_dict.items():
        if isinstance(repeat_division_num, dict):
            this_result_repeat_division = repeat_division_num[result_label]
        elif isinstance(repeat_division_num, int):
            this_result_repeat_division = repeat_division_num
        elif repeat_division_num is None:
            this_result_repeat_division = None
        else:
            raise TypeError
        if isinstance(each_optimization_num, dict):
            this_result_optimization_num = each_optimization_num[result_label]
        elif isinstance(each_optimization_num, int):
            this_result_optimization_num = each_optimization_num
        elif each_optimization_num is None:
            this_result_optimization_num = None
        else:
            raise TypeError
        solver_obj = solver_dict[result_label]
        solution_array = final_solution_data_dict[result_label]
        predicted_mid_data_dict = final_predicted_mid_data_dict[result_label]
        new_loss_list = []
        new_solution_list = []
        new_mid_data_dict = {}
        total_solution_num = len(loss_array)
        for index_start, index_end in index_start_end_generator(
                total_solution_num, this_result_repeat_division, this_result_optimization_num):
            this_repeat_loss_array = loss_array[index_start:index_end]
            filtered_index = select_solutions(this_repeat_loss_array, loss_percentile, select_num, index_start)
            mean_solution_array = np.mean(solution_array[filtered_index, :], 0)
            calculated_loss = solver_obj.obj_eval(mean_solution_array)
            calculated_mid_data_dict = solver_obj.predict(mean_solution_array)
            new_solution_list.append(mean_solution_array)
            new_loss_list.append(calculated_loss)
            for mid_title, predicted_mid_array in calculated_mid_data_dict.items():
                processed_mid_title = mid_name_process(mid_title)
                if processed_mid_title not in new_mid_data_dict:
                    new_mid_data_dict[processed_mid_title] = []
                new_mid_data_dict[processed_mid_title].append(predicted_mid_array)
        new_loss_data_dict[result_label] = np.array(new_loss_list)
        new_final_solution_data_dict[result_label] = np.array(new_solution_list)
        new_final_predicted_mid_data_dict[result_label] = new_mid_data_dict

    return new_final_solution_data_dict, new_loss_data_dict, new_final_predicted_mid_data_dict


def experimental_mid_prediction_plotting_func_template(
        tissue_mid_name_list_dict_constructor, result_process_information_dict_analysis,
        final_information_dict, mid_name_tissue_name_dict, solver_dict,
        legend_color_dict=None, **figure_param_dict):
    """
    Template function for experimental MID prediction plotting.
    Returns a function that performs the actual plotting.

    Similar pattern to experimental_data_plotting_func_template.
    """
    default_experimental_transparency = 0.3

    def preprocess_predicted_data(
            complex_predicted_data_dict, final_target_experimental_mid_data_dict, subset_index_dict):
        """Preprocess all result label data for MID prediction plotting."""
        tmp_formatted_mean_data_dict = {}
        tmp_formatted_std_data_dict = {}
        data_len_dict = {}
        subplot_name_dict = {}
        tmp_color_dict = {}

        tissue_mid_name_list_dict, display_mid_name_dict, display_tissue_name_dict = tissue_mid_name_list_dict_constructor(
            result_process=True)

        for (
                result_label, major_key, minor_key_list, minor_key_str, current_color, order_index
        ) in result_process_information_dict_analysis(final_information_dict, target=Keywords.mid_prediction):
            # Filter out information for results that are missing in the predicted data
            # (e.g. failed runs with 0 solutions)
            if result_label not in complex_predicted_data_dict:
                continue
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
            experimental_key_str = f'{minor_key_str}_experimental'
            if major_key not in tmp_color_dict:
                tmp_color_dict[major_key] = {}
            tmp_color_dict[major_key].update({
                predicted_key_str: (order_index, current_color),
                experimental_key_str: (order_index + 0.5, current_color.transparency_mix(default_experimental_transparency)),
            })
            current_target_experimental_mid_data_dict = final_target_experimental_mid_data_dict[result_label]
            # TODO: Need to mark the existing but manually excluded MFA data, similar with experimental_data_plotting.

            for mid_title, current_predicted_data_array_list in result_specific_predicted_data_dict.items():
                if mid_title in mid_name_tissue_name_dict:
                    tissue_name, raw_mid_name = mid_name_tissue_name_dict[mid_title]
                else:
                    from ..core.common.functions import tissue_name_breakdown
                    tissue_name, raw_mid_name = tissue_name_breakdown(mid_title)
                    mid_name_tissue_name_dict[mid_title] = (tissue_name, raw_mid_name)
                current_experimental_name = current_emu_name_experimental_name_dict[mid_title]
                updated_tissue_specific_mid_name = tissue_specific_name_constructor(current_experimental_name, tissue_name)
                current_experimental_mid = current_solver_target_experimental_mid_dict[mid_title]
                if major_key not in tmp_formatted_mean_data_dict:
                    tmp_formatted_mean_data_dict[major_key] = {}
                    tmp_formatted_std_data_dict[major_key] = {}
                current_tmp_formatted_mean_data_dict = tmp_formatted_mean_data_dict[major_key]
                current_tmp_formatted_std_data_dict = tmp_formatted_std_data_dict[major_key]
                if updated_tissue_specific_mid_name not in current_tmp_formatted_mean_data_dict:
                    current_tmp_formatted_mean_data_dict[updated_tissue_specific_mid_name] = {}
                    current_tmp_formatted_std_data_dict[updated_tissue_specific_mid_name] = {}
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
                current_tmp_formatted_mean_data_dict[updated_tissue_specific_mid_name][predicted_key_str] = (
                    order_index, target_predicted_data_array.mean(axis=0))
                current_tmp_formatted_std_data_dict[updated_tissue_specific_mid_name][predicted_key_str] = (
                    order_index, target_predicted_data_array.std(axis=0))
                current_tmp_formatted_mean_data_dict[updated_tissue_specific_mid_name][experimental_key_str] = (
                    order_index + 0.5, current_experimental_mid
                )

        formatted_mean_data_dict = convert_tmp_data_to_final_data_dict(tmp_formatted_mean_data_dict)
        formatted_std_data_dict = convert_tmp_data_to_final_data_dict(tmp_formatted_std_data_dict)
        complete_color_dict = {}
        for major_key, each_major_key_tmp_color_dict in tmp_color_dict.items():
            complete_color_dict[major_key] = sort_dict_and_remove_sorting_index(each_major_key_tmp_color_dict)
        return (
            formatted_mean_data_dict, formatted_std_data_dict, data_len_dict, subplot_name_dict, complete_color_dict,
            tissue_mid_name_list_dict, display_mid_name_dict, display_tissue_name_dict
        )

    def mid_prediction_plotting(
            experiment_name, complex_predicted_data_dict, final_target_experimental_mid_data_dict,
            prediction_output_direct, subset_index_dict=None):
        """Main plotting function that processes data and generates figures."""
        (
            formatted_mean_data_dict, formatted_std_data_dict, data_len_dict, subplot_name_dict, complete_color_dict,
            tissue_mid_name_list_dict, display_mid_name_dict, display_tissue_name_dict
        ) = preprocess_predicted_data(
            complex_predicted_data_dict, final_target_experimental_mid_data_dict, subset_index_dict)

        mean_data_dict = {}
        std_data_dict = {}
        mid_name_list_dict = {}
        current_legend_color_dict = legend_color_dict

        for major_key, current_formatted_mean_data in formatted_mean_data_dict.items():
            current_formatted_std_data = formatted_std_data_dict[major_key]
            # Import helper function for tissue-specific metabolite MID list construction
            current_mid_name_list, current_mean_data_dict, current_std_data_dict = tissue_specific_metabolite_mid_list_constructor(
                current_formatted_mean_data, current_formatted_std_data, data_len_dict, tissue_mid_name_list_dict)
            mean_data_dict[major_key] = current_formatted_mean_data
            std_data_dict[major_key] = current_formatted_std_data
            mid_name_list_dict[major_key] = current_mid_name_list

            target_row_num = len(current_mid_name_list)
            target_col_num = max([len(each_row_mid_name_list) for each_row_mid_name_list in current_mid_name_list])
            if current_legend_color_dict is not None:
                major_legend_color_dict = current_legend_color_dict[major_key]
            else:
                major_legend_color_dict = None
            color_dict = complete_color_dict[major_key]

            figure_plotting.mid_comparison_plotting(
                figure_output_direct=prediction_output_direct, figure_name=f'predicted_mid_comparison_{major_key}',
                complete_data_dict=current_mean_data_dict, target_emu_name_nested_list=current_mid_name_list,
                target_row_num=target_row_num, target_col_num=target_col_num, error_bar_data_dict=current_std_data_dict,
                color_dict=color_dict, legend_color_dict=major_legend_color_dict, y_lim=(0, 1),
                **figure_param_dict)

        figure_raw_data = FigureData(FigureDataKeywords.mid_comparison, experiment_name)
        figure_raw_data.save_data(
            mean_data_dict=mean_data_dict,
            std_data_dict=std_data_dict,
            mid_name_list=mid_name_list_dict,
            subplot_name_dict=subplot_name_dict,
        )

    return mid_prediction_plotting


def common_flux_comparison_func(
        current_important_flux_list, common_flux_name_index_dict, current_data_array, data_dict_for_plotting, key_name,
        order_index, reversible_flux_title_constructor=None, extra_group_key_name=None):
    def default_reversible_flux_title_constructor(flux_name_0, flux_name_1):
        return f'{flux_name_0} - {flux_name_1}'

    if reversible_flux_title_constructor is None:
        reversible_flux_title_constructor = default_reversible_flux_title_constructor
    for flux_name in current_important_flux_list:
        if isinstance(flux_name, str):
            flux_title = flux_name
            flux_index = common_flux_name_index_dict[flux_name]
            calculated_flux_array = current_data_array[:, flux_index]
        elif isinstance(flux_name, tuple | list):
            if callable(flux_name[1]):
                flux_title, flux_func = flux_name
                flux_name_value_dict = {
                    tmp_flux_name: current_data_array[:, flux_index]
                    for tmp_flux_name, flux_index in common_flux_name_index_dict.items()}
                calculated_flux_array = flux_func(flux_name_value_dict)
            else:
                # flux_title = '{} - {}'.format(flux_name[0], flux_name[1])
                flux_title = reversible_flux_title_constructor(*flux_name)
                flux_index1 = common_flux_name_index_dict[flux_name[0]]
                flux_index2 = common_flux_name_index_dict[flux_name[1]]
                calculated_flux_array = (
                    current_data_array[:, flux_index1] - current_data_array[:, flux_index2])
        else:
            raise ValueError()
        if flux_title not in data_dict_for_plotting:
            data_dict_for_plotting[flux_title] = {}
        if extra_group_key_name is not None:
            if extra_group_key_name not in data_dict_for_plotting[flux_title]:
                data_dict_for_plotting[flux_title][extra_group_key_name] = {}
            data_dict_for_plotting[flux_title][extra_group_key_name][key_name] = (order_index, calculated_flux_array)
        else:
            data_dict_for_plotting[flux_title][key_name] = (order_index, calculated_flux_array)


def obtain_current_specific_flux_array(
        buffer_dict, final_flux_name_index_dict, final_solution_data_dict, subset_index_dict, result_label,
        *tissue_flux_name_list):
    def filter_subset_index(
            _raw_solution_array, _subset_index_dict, _result_label, _flux_index):
        if _subset_index_dict is not None:
            _solution_array = _raw_solution_array[_subset_index_dict[_result_label], _flux_index]
        else:
            _solution_array = _raw_solution_array[:, _flux_index]
        return _solution_array

    def select_flux(
            _current_tissue_flux_name_index_dict, _current_raw_solution_array, _subset_index_dict,
            _result_label, _tissue_flux_name):
        if _tissue_flux_name is None:
            return None
        elif (_result_label, _tissue_flux_name) in buffer_dict:
            return buffer_dict[(_result_label, _tissue_flux_name)]
        else:
            try:
                _current_flux_index = _current_tissue_flux_name_index_dict[_tissue_flux_name]
            except KeyError:
                _current_flux_array = None
            else:
                _current_flux_array = filter_subset_index(
                    _current_raw_solution_array, _subset_index_dict, _result_label, _current_flux_index)
            buffer_dict[(_result_label, _tissue_flux_name)] = _current_flux_array
            return _current_flux_array

    def obtain_one_result_label(
            _final_flux_name_index_dict, _final_solution_data_dict, _subset_index_dict,
            _result_label, _tissue_flux_name_list):
        _current_tissue_flux_name_index_dict = _final_flux_name_index_dict[_result_label]
        _current_raw_solution_array = _final_solution_data_dict[_result_label]
        if len(_tissue_flux_name_list) == 0:
            _tissue_flux_name_list = _current_tissue_flux_name_index_dict.keys()
        _flux_array_list = []
        for _tissue_flux_name in _tissue_flux_name_list:
            _current_flux_array = select_flux(
                _current_tissue_flux_name_index_dict, _current_raw_solution_array, _subset_index_dict,
                _result_label, _tissue_flux_name
            )
            _flux_array_list.append(_current_flux_array)
        return _flux_array_list, _tissue_flux_name_list


    if isinstance(result_label, str):
        flux_array_list, new_tissue_flux_name_list = obtain_one_result_label(
            final_flux_name_index_dict, final_solution_data_dict, subset_index_dict, result_label,
            tissue_flux_name_list)
    elif isinstance(result_label, tuple | list):
        total_num_list = []
        tissue_flux_num = len(tissue_flux_name_list)
        all_flux_array_list = []
        new_tissue_flux_name_list = None
        for each_result_label in result_label:
            this_result_label_flux_array_list, this_result_label_new_tissue_flux_name_list = obtain_one_result_label(
                final_flux_name_index_dict, final_solution_data_dict, subset_index_dict, each_result_label,
                tissue_flux_name_list)
            if tissue_flux_num == 0:
                tissue_flux_num = len(this_result_label_new_tissue_flux_name_list)
                new_tissue_flux_name_list = this_result_label_new_tissue_flux_name_list
            else:
                assert tissue_flux_num == len(this_result_label_new_tissue_flux_name_list)
            if len(all_flux_array_list) == 0:
                all_flux_array_list = [[] for _ in range(tissue_flux_num)]
                total_num_list = [0 for _ in range(tissue_flux_num)]
            for flux_order_index in range(tissue_flux_num):
                each_flux_array = this_result_label_flux_array_list[flux_order_index]
                this_flux_each_flux_array_list = all_flux_array_list[flux_order_index]
                if each_flux_array is not None:
                    this_flux_each_flux_array_list.append(each_flux_array)
                    this_solution_num = each_flux_array.shape[0]
                    total_num_list[flux_order_index] += this_solution_num
        flux_array_list = [
            np.concatenate(each_flux_array_list, axis=0) if len(each_flux_array_list) > 0 else None
            for each_flux_array_list in all_flux_array_list
        ]
    else:
        raise ValueError()

    return flux_array_list, new_tissue_flux_name_list


def process_each_result_label(
        buffer_dict, result_label, final_solution_data_dict, final_flux_name_index_dict, subset_index_dict,
        all_tissue_flux_name_list_dict, display_flux_dict, display_tissue_dict,
        final_flux_name_nested_list, tissue_specific_display_flux_dict):
    def tissue_flux_label_generator(
            _tissue_name, _flux_name_without_tissue_0, _flux_name_without_tissue_1=None):
        if _flux_name_without_tissue_1 is None:
            _tissue_flux_name = tissue_specific_name_constructor(_flux_name_without_tissue_0, _tissue_name)
            _flux_label = _tissue_flux_name
        else:
            _flux_label_without_tissue = f'{_flux_name_without_tissue_0} - {_flux_name_without_tissue_1}'
            _flux_label = tissue_specific_name_constructor(_flux_label_without_tissue, _tissue_name)
            _tissue_flux_name_0 = tissue_specific_name_constructor(_flux_name_without_tissue_0, _tissue_name)
            _tissue_flux_name_1 = tissue_specific_name_constructor(_flux_name_without_tissue_1, _tissue_name)
            _tissue_flux_name = (_tissue_flux_name_0, _tissue_flux_name_1)
        return _flux_label, _tissue_flux_name

    def tissue_display_flux_name_generator(
            _display_flux_dict, _display_tissue_dict,
            _tissue_name, _flux_name_without_tissue_0, _flux_name_without_tissue_1=None,):
        if _tissue_name in _display_tissue_dict:
            _display_tissue = _display_tissue_dict[_tissue_name]
        else:
            _display_tissue = _tissue_name
        if _flux_name_without_tissue_1 is None:
            if _flux_name_without_tissue_0 in _display_flux_dict:
                _display_flux_name = _display_flux_dict[_flux_name_without_tissue_0]
            else:
                _display_flux_name = _flux_name_without_tissue_0
        else:
            flux_pair_tuple = (_flux_name_without_tissue_0, _flux_name_without_tissue_1)
            if flux_pair_tuple in _display_flux_dict:
                _display_flux_name = _display_flux_dict[flux_pair_tuple]
            elif _flux_name_without_tissue_0 in _display_flux_dict:
                _display_flux_name = f'{_display_flux_dict[_flux_name_without_tissue_0]}_net'
            else:
                _display_flux_name = ' - '.join(flux_pair_tuple)
        _display_tissue_flux_name = f'{_display_tissue}: {_display_flux_name}'
        return _display_tissue_flux_name

    def reversible_flux_label_constructor(_flux_name_tuple, _tissue_name):
        if isinstance(_flux_name_tuple[1], str):
            # flux_name_tuple = (flux_name_without_tissue_0, flux_name_without_tissue_1)
            _flux_label, (_tissue_flux_name_0, _tissue_flux_name_1) = tissue_flux_label_generator(
                _tissue_name, *_flux_name_tuple)
            _display_tissue_flux_name = tissue_display_flux_name_generator(
                display_flux_dict, display_tissue_dict,
                _tissue_name, *_flux_name_tuple)
        elif isinstance(_flux_name_tuple[1], tuple | list):
            assert isinstance(_flux_name_tuple[0], str)
            assert len(_flux_name_tuple[1]) == 2
            (_flux_label, _flux_pair_tuple) = _flux_name_tuple
            _, (_tissue_flux_name_0, _tissue_flux_name_1) = tissue_flux_label_generator(
                _tissue_name, *_flux_pair_tuple)
            _display_tissue_flux_name = tissue_display_flux_name_generator(
                display_flux_dict, display_tissue_dict,
                _tissue_name, *_flux_pair_tuple)
        else:
            raise ValueError()
        return _flux_label, _display_tissue_flux_name, _tissue_flux_name_0, _tissue_flux_name_1

    result_flux_dict = {}
    if len(final_flux_name_nested_list) == 0:
        add_flux_name_state = True
    else:
        add_flux_name_state = False

    total_row_index = 0
    for tissue_name, each_tissue_flux_name_nested_list in all_tissue_flux_name_list_dict.items():
        for each_tissue_flux_name_list in each_tissue_flux_name_nested_list:
            each_row_flux_label_list = []
            for flux_name_without_tissue in each_tissue_flux_name_list:
                if isinstance(flux_name_without_tissue, str):
                    flux_label, tissue_flux_name = tissue_flux_label_generator(tissue_name, flux_name_without_tissue)
                    (calculated_flux_array, ), _ = obtain_current_specific_flux_array(
                        buffer_dict, final_flux_name_index_dict, final_solution_data_dict,
                        subset_index_dict, result_label, tissue_flux_name)
                    if calculated_flux_array is None:
                        flux_label = None
                        flux_display_name = None
                    else:
                        flux_display_name = tissue_display_flux_name_generator(
                            display_flux_dict, display_tissue_dict, tissue_name, flux_name_without_tissue)
                elif isinstance(flux_name_without_tissue, tuple | list):
                    if len(flux_name_without_tissue) == 2:
                        # flux_title = '{} - {}'.format(flux_name[0], flux_name[1])
                        (
                            flux_label, flux_display_name, tissue_flux_name_0, tissue_flux_name_1
                        ) = reversible_flux_label_constructor(flux_name_without_tissue, tissue_name)

                        (calculated_flux_array_0, calculated_flux_array_1), _ = obtain_current_specific_flux_array(
                            buffer_dict, final_flux_name_index_dict, final_solution_data_dict,
                            subset_index_dict, result_label, tissue_flux_name_0, tissue_flux_name_1)
                        if calculated_flux_array_0 is not None and calculated_flux_array_1 is not None:
                            calculated_flux_array = calculated_flux_array_0 - calculated_flux_array_1
                        else:
                            flux_label = None
                            flux_display_name = None
                            calculated_flux_array = None
                    elif len(flux_name_without_tissue) == 3 and callable(flux_name_without_tissue[2]):
                        flux_label, flux_display_name, flux_func = flux_name_without_tissue
                        flux_array_list, new_tissue_flux_name_list = obtain_current_specific_flux_array(
                            buffer_dict, final_flux_name_index_dict, final_solution_data_dict,
                            subset_index_dict, result_label, )
                        tissue_flux_name_value_dict = {
                            tmp_tissue_flux_name: flux_array
                            for tmp_tissue_flux_name, flux_array in zip(new_tissue_flux_name_list, flux_array_list)
                        }
                        calculated_flux_array = flux_func(tissue_flux_name_value_dict)
                    else:
                        raise ValueError()
                else:
                    raise ValueError()
                if flux_label is not None:
                    assert flux_label not in result_flux_dict
                    result_flux_dict[flux_label] = calculated_flux_array
                    if add_flux_name_state:
                        tissue_specific_display_flux_dict[flux_label] = flux_display_name
                    else:
                        assert tissue_specific_display_flux_dict[flux_label] == flux_display_name
                each_row_flux_label_list.append(flux_label)
            if add_flux_name_state:
                final_flux_name_nested_list.append(each_row_flux_label_list)
            else:
                existing_row_flux_label_list = final_flux_name_nested_list[total_row_index]
                assert len(existing_row_flux_label_list) == len(each_row_flux_label_list)
                new_row_flux_label_list = []
                for col_index, (existing_flux_label, new_flux_label) in enumerate(
                        zip(existing_row_flux_label_list, each_row_flux_label_list)):
                    if existing_flux_label is not None and new_flux_label is not None:
                        assert existing_flux_label == new_flux_label
                    elif existing_flux_label is None and new_flux_label is not None:
                        new_list_len = len(new_row_flux_label_list)
                        if new_list_len < col_index:
                            new_row_flux_label_list.extend(existing_row_flux_label_list[new_list_len:col_index])
                            new_row_flux_label_list.append(new_flux_label)
                if len(new_row_flux_label_list) > 0:
                    new_list_len = len(new_row_flux_label_list)
                    new_row_flux_label_list.extend(existing_row_flux_label_list[new_list_len:])
                    final_flux_name_nested_list[total_row_index] = new_row_flux_label_list
            total_row_index += 1

    return (result_flux_dict, )


def flux_comparison_plotting_func_template(
        tissue_flux_name_list_dict_constructor, flux_name_tissue_name_dict, result_process_information_dict_analysis,
        reversible_flux_pair_dict_generator, tissue_specific_reversible_flux_title_constructor_generator,
        **figure_param_dict):
    """
    Template function for flux comparison plotting.
    Returns a function that performs the actual plotting.

    Similar pattern to experimental_data_plotting_func_template.
    """

    def flux_comparison_plotting(
            experiment_name, solution_array_dict, final_flux_name_index_dict, final_information_dict,
            flux_comparison_output_direct, subset_index_dict=None):
        """Main plotting function that processes data and generates figures."""
        # preprocessed_flux_data_dict, all_tissue_flux_name_list_dict, display_flux_dict = preprocess_flux_data(
        #     solution_array_dict, final_flux_name_index_dict, subset_index_dict)

        tissue_tissue_flux_name_list_dict_info = tissue_flux_name_list_dict_constructor()
        if len(tissue_tissue_flux_name_list_dict_info) == 5:
            (
                all_tissue_flux_name_list_dict, display_flux_dict, display_tissue_dict,
                specific_tissue_reversible_flux_pair_dict, flux_process_map_func
            ) = tissue_tissue_flux_name_list_dict_info
        elif len(tissue_tissue_flux_name_list_dict_info) == 4:
            (
                all_tissue_flux_name_list_dict, display_flux_dict, display_tissue_dict,
                specific_tissue_reversible_flux_pair_dict,
            ) = tissue_tissue_flux_name_list_dict_info
            flux_process_map_func = None
        elif len(tissue_tissue_flux_name_list_dict_info) == 3:
            (
                all_tissue_flux_name_list_dict, display_flux_dict,
                specific_tissue_reversible_flux_pair_dict,
            ) = tissue_tissue_flux_name_list_dict_info
            display_tissue_dict = {}
            flux_process_map_func = None
        else:
            raise ValueError()

        final_flux_name_nested_list = []
        tissue_specific_display_flux_dict = {}
        processed_result_dict = {}
        extra_group_key_name = ''
        tmp_organized_mouse_data_dict = {}
        tmp_color_dict = {}
        tmp_name_dict = {}
        buffer_dict = {}
        for (
                result_label, major_key, minor_key_list, minor_key_str, current_color, order_index
        ) in result_process_information_dict_analysis(final_information_dict, target=Keywords.flux_comparison):
            if result_label in processed_result_dict:
                processed_result = processed_result_dict[result_label]
            else:
                processed_result = process_each_result_label(
                    buffer_dict, result_label, solution_array_dict, final_flux_name_index_dict, subset_index_dict,
                    all_tissue_flux_name_list_dict, display_flux_dict, display_tissue_dict,
                    final_flux_name_nested_list, tissue_specific_display_flux_dict)
                processed_result_dict[result_label] = processed_result
            (processed_flux_value_dict, ) = processed_result

            if major_key not in tmp_organized_mouse_data_dict:
                tmp_organized_mouse_data_dict[major_key] = {}
                tmp_color_dict[major_key] = {}
                tmp_name_dict[major_key] = {}
            tmp_color_dict[major_key][minor_key_str] = (order_index, current_color)
            tmp_name_dict[major_key][minor_key_str] = (order_index, minor_key_str)
            for flux_title, each_flux_data_array in processed_flux_value_dict.items():
                if flux_process_map_func is not None:
                    mapped_each_flux_data_array = flux_process_map_func(each_flux_data_array, major_key, minor_key_str)
                else:
                    mapped_each_flux_data_array = each_flux_data_array
                current_flux_tuple = (order_index, mapped_each_flux_data_array)
                if flux_title not in tmp_organized_mouse_data_dict[major_key]:
                    tmp_organized_mouse_data_dict[major_key][flux_title] = {}
                current_parent_data_dict = tmp_organized_mouse_data_dict[major_key][flux_title]
                if extra_group_key_name not in current_parent_data_dict:
                    current_parent_data_dict[extra_group_key_name] = {}
                current_parent_data_dict[extra_group_key_name][minor_key_str] = current_flux_tuple

        organized_mouse_data_dict = {}
        for major_key, each_major_key_tmp_data_dict in tmp_organized_mouse_data_dict.items():
            ordered_color_dict = sort_dict_and_remove_sorting_index(tmp_color_dict[major_key])
            ordered_name_dict = sort_dict_and_remove_sorting_index(tmp_name_dict[major_key])
            current_organized_mouse_data_dict = {}
            for flux_name, each_flux_data_dict in each_major_key_tmp_data_dict.items():
                current_organized_mouse_data_dict[flux_name] = {
                    extra_group_key_name: sort_dict_and_remove_sorting_index(
                        each_flux_data_dict[extra_group_key_name])}
            figure_plotting.distribution_comparison_plotting(
                figure_output_direct=flux_comparison_output_direct, figure_name=f'flux_comparison_{major_key}',
                complete_data_dict=current_organized_mouse_data_dict, flux_name_nested_list=final_flux_name_nested_list,
                display_flux_dict=tissue_specific_display_flux_dict, name_dict=ordered_name_dict,
                color_dict=ordered_color_dict, **figure_param_dict)
            organized_mouse_data_dict[major_key] = current_organized_mouse_data_dict

        figure_raw_data = FigureData(FigureDataKeywords.flux_comparison, experiment_name)
        figure_raw_data.save_data(
            organized_diet_mouse_data_dict=organized_mouse_data_dict,)
        return organized_mouse_data_dict

    return flux_comparison_plotting

