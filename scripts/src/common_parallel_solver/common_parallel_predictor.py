from .packages import np, mp, it, threadpool_limits
from .config import (
    Keywords, random_seed, specific_solver_constructor, base_solver_constructor,
    parameter_extract, split_total_num_to_process, ProgressBarExtraProcess)

from ..core.common.functions import (
    compartmental_mid_name_constructor, mid_name_process)
from ..core.solver.slsqp_numba_solver.solver_class import SLSQPNumbaSolver


def parallel_parameter_generator(
        flux_vector_array_list, base_solver_obj_list, mfa_solver_option_dict_list, flux_vector_index_num_nested_list,
        thread_num_constraint, predict_all_mid, send_pipe):

    for solver_index, (
            flux_vector_array, base_solver_obj, mfa_solver_option_dict, flux_vector_index_num_list
    ) in enumerate(zip(
        flux_vector_array_list, base_solver_obj_list, mfa_solver_option_dict_list, flux_vector_index_num_nested_list
    )):
        flux_vector_start_index = 0
        for flux_vector_num in flux_vector_index_num_list:
            current_flux_vector_array = flux_vector_array[flux_vector_start_index:flux_vector_start_index + flux_vector_num]
            parameter_list = (
                solver_index, base_solver_obj, mfa_solver_option_dict, current_flux_vector_array, predict_all_mid,
                thread_num_constraint, send_pipe)
            flux_vector_start_index += flux_vector_num
            yield parameter_list


def pure_predictor(slsqp_solver_obj, flux_vector, predict_all_mid=False):
    raw_predicted_data_dict = slsqp_solver_obj.predict(flux_vector)
    if predict_all_mid:
        raw_all_predicted_data_dict = slsqp_solver_obj.predict_all_target(flux_vector)
    else:
        raw_all_predicted_data_dict = {}
    return raw_predicted_data_dict, raw_all_predicted_data_dict


def parallel_single_predictor(parameter_list):
    (
        solver_index, base_solver_obj, mfa_solver_option_dict, flux_vector_array, predict_all_mid,
        thread_num_constraint, send_pipe) = parameter_list
    slsqp_solver_obj = SLSQPNumbaSolver(base_solver_obj, solver_option_dict=mfa_solver_option_dict)
    slsqp_solver_obj.initialize_solver()
    raw_predicted_data_dict_list = []
    raw_all_predicted_data_dict_list = []
    for flux_vector in flux_vector_array:
        with threadpool_limits(limits=thread_num_constraint):
            raw_predicted_data_dict, raw_all_predicted_data_dict = pure_predictor(
                slsqp_solver_obj, flux_vector, predict_all_mid)
        raw_predicted_data_dict_list.append(raw_predicted_data_dict)
        if predict_all_mid:
            raw_all_predicted_data_dict_list.append(raw_all_predicted_data_dict)
        if send_pipe is not None:
            send_pipe.send(1)
    return solver_index, raw_predicted_data_dict_list, raw_all_predicted_data_dict_list


def common_parallel_predictor(
        parameter_list_iter, total_solver_num, processes_num=6, parallel_test=False, predict_all_mid=False,
        **other_parameters):
    def process_result(current_raw_result):
        (
            solver_index, raw_predicted_data_dict_list, raw_all_predicted_data_dict_list) = current_raw_result
        folded_raw_predicted_data_dict_nested_list[solver_index].append(raw_predicted_data_dict_list)
        if predict_all_mid:
            folded_raw_all_predicted_data_dict_nested_list[solver_index].append(raw_all_predicted_data_dict_list)

    folded_raw_predicted_data_dict_nested_list = [[] for _ in range(total_solver_num)]
    folded_raw_all_predicted_data_dict_nested_list = [[] for _ in range(total_solver_num)]

    if parallel_test:
        for parameter_list in parameter_list_iter:
            raw_result = parallel_single_predictor(parameter_list)
            process_result(raw_result)
    else:
        with mp.Pool(processes=processes_num) as pool:
            raw_result_iter = pool.imap(parallel_single_predictor, parameter_list_iter)
            for raw_result in raw_result_iter:
                process_result(raw_result)

    all_final_raw_predicted_mid_data_dict_list = []
    for folded_raw_predicted_data_dict_list in folded_raw_predicted_data_dict_nested_list:
        raw_long_name_predicted_mid_data_dict = {}
        for raw_predicted_data_dict_list in folded_raw_predicted_data_dict_list:
            for raw_predicted_data_dict in raw_predicted_data_dict_list:
                for raw_mid_name, mid_data_array in raw_predicted_data_dict.items():
                    if raw_mid_name not in raw_long_name_predicted_mid_data_dict:
                        raw_long_name_predicted_mid_data_dict[raw_mid_name] = []
                    raw_long_name_predicted_mid_data_dict[raw_mid_name].append(mid_data_array)
        final_raw_predicted_mid_data_dict = {
            mid_name_process(raw_mid_name): np.array(mid_data_array_list)
            for raw_mid_name, mid_data_array_list in raw_long_name_predicted_mid_data_dict.items()}
        all_final_raw_predicted_mid_data_dict_list.append(final_raw_predicted_mid_data_dict)

    all_final_raw_all_predicted_mid_data_dict_list = []
    for folded_raw_all_predicted_data_dict_list in folded_raw_all_predicted_data_dict_nested_list:
        raw_list_all_predicted_mid_data_dict = {}
        for raw_all_predicted_data_dict_list in folded_raw_all_predicted_data_dict_list:
            for raw_all_predicted_data_dict in raw_all_predicted_data_dict_list:
                for mid_name, mid_data_array in raw_all_predicted_data_dict.items():
                    if mid_name not in raw_list_all_predicted_mid_data_dict:
                        raw_list_all_predicted_mid_data_dict[mid_name] = []
                    raw_list_all_predicted_mid_data_dict[mid_name].append(mid_data_array)
        final_raw_all_predicted_mid_data_dict = {
            raw_mid_name: np.array(mid_data_array_list)
            for raw_mid_name, mid_data_array_list in raw_list_all_predicted_mid_data_dict.items()}
        all_final_raw_all_predicted_mid_data_dict_list.append(final_raw_all_predicted_mid_data_dict)

    return all_final_raw_predicted_mid_data_dict_list, all_final_raw_all_predicted_mid_data_dict_list


def base_solver_obj_generator(
        flux_vector_array_list, mfa_tuple_list, slsqp_solver_obj_list, processes_num, parallel_test):
    total_solver_num = len(flux_vector_array_list)
    base_solver_obj_for_parallel_list = []
    mfa_solver_option_dict_list = []
    if slsqp_solver_obj_list is not None:
        assert len(slsqp_solver_obj_list) == total_solver_num
        for slsqp_solver_obj in slsqp_solver_obj_list:
            slsqp_solver_option_dict = slsqp_solver_obj.solver_option_dict
            base_solver_obj_for_parallel = slsqp_solver_obj.__copy__()
            base_solver_obj_for_parallel.base_initialize_solver()
            base_solver_obj_for_parallel_list.append(base_solver_obj_for_parallel)
            mfa_solver_option_dict_list.append(slsqp_solver_option_dict)
    else:
        assert mfa_tuple_list is not None and len(mfa_tuple_list) == total_solver_num
        for mfa_model, mfa_data, mfa_config in mfa_tuple_list:
            base_solver_obj_for_parallel = base_solver_constructor(mfa_model, mfa_data, mfa_config)
            base_solver_obj_for_parallel.base_initialize_solver()
            base_solver_obj_for_parallel_list.append(base_solver_obj_for_parallel)
            mfa_solver_option_dict_list.append(mfa_config.solver_config_dict)
    flux_vector_index_num_nested_list = []
    total_flux_size = 0
    for flux_vector_array in flux_vector_array_list:
        current_flux_size = flux_vector_array.shape[0]
        total_flux_size += current_flux_size
        current_flux_vector_index_num_list = split_total_num_to_process(
            current_flux_size, processes_num)
        flux_vector_index_num_nested_list.append(current_flux_vector_index_num_list)
    updated_parallel_test = parallel_test or total_flux_size < 5

    return (
        base_solver_obj_for_parallel_list, mfa_solver_option_dict_list, total_solver_num, total_flux_size,
        updated_parallel_test, flux_vector_index_num_nested_list)


def multi_solver_parallel_predictor(
        flux_vector_array_list, mfa_tuple_list=None, slsqp_solver_obj_list=None,
        processes_num=6, predict_all_mid=False, parallel_test=False, display_progress_bar=False, name=None):
    thread_num_constraint = max(4, int(24 / processes_num))
    (
        base_solver_obj_for_parallel_list, mfa_solver_option_dict_list, total_solver_num, total_size,
        updated_parallel_test, flux_vector_index_num_nested_list,
     ) = base_solver_obj_generator(
        flux_vector_array_list, mfa_tuple_list, slsqp_solver_obj_list, processes_num, parallel_test)
    if name is not None:
        display_title = f'MID prediction for {name}'
    else:
        display_title = 'MID prediction'
    with ProgressBarExtraProcess(
            display_progress_bar=display_progress_bar, target_size=total_size, display_title=display_title,
    ) as progress_bar:
        parameter_list_iter = parallel_parameter_generator(
            flux_vector_array_list, base_solver_obj_for_parallel_list, mfa_solver_option_dict_list,
            flux_vector_index_num_nested_list, thread_num_constraint=thread_num_constraint,
            predict_all_mid=predict_all_mid, send_pipe=progress_bar.send_pipe)
        (
            all_final_raw_predicted_mid_data_dict_list, all_final_raw_all_predicted_mid_data_dict_list
        ) = common_parallel_predictor(
            parameter_list_iter, total_solver_num=total_solver_num, processes_num=processes_num,
            parallel_test=updated_parallel_test, predict_all_mid=predict_all_mid,)
    return all_final_raw_predicted_mid_data_dict_list, all_final_raw_all_predicted_mid_data_dict_list

