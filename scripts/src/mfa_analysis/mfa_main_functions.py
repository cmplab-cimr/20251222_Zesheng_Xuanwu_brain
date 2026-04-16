import sys

from ..core.solver.solver_construction_functions.solver_constructor import common_solver_constructor
from ..common_parallel_solver.common_parallel_solver import common_solver

from ..io_functions.result_output_functions import solver_output

from .config import (
    ParallelSolverKeywords, ExperimentName, MFARunningMode, Direct, default_data_model_name)
from .result_processing_and_plotting import (
    FinalResult, OutOfOrderFinalResult, experimental_mid_and_raw_data_plotting)


def solver_dict_constructor(parameter_label_content_dict):
    target_solver_dict = {}
    same_model_dict = {}
    same_data_dict = {}
    previous_label = None
    for result_label, (_, (mfa_model, mfa_data, mfa_config), _, _) in parameter_label_content_dict.items():
        current_solver_obj = common_solver_constructor(
            mfa_model, mfa_data, mfa_config, name=result_label, verbose=False)
        target_solver_dict[result_label] = current_solver_obj
        if previous_label is None:
            same_model_dict[result_label] = False
            previous_label = result_label
        else:
            same_model_dict[result_label] = True
        if result_label.startswith(previous_label) and result_label.endswith(ParallelSolverKeywords.unoptimized):
            same_data_dict[result_label] = True
        else:
            same_data_dict[result_label] = False
    return target_solver_dict, same_model_dict, same_data_dict


def old_model_data_config_generator(
        data_model_object, user_defined_mfa_settings, running_mode, experiment_name, dataset_id,
        test_mode, docker_mode, suffix, parallel_num):
    # 用 mfa_data_dict 的 keys 来构建模型映射，且所有 result_label 使用同一模型
    (
        mfa_data_dict, experimental_result_information_dict_analysis,
        experimental_tissue_mid_name_list_dict_constructor,
        result_information_dict, experimental_figure_config_dict,
    ) = data_model_object.return_data_content(
        running_mode=running_mode, experiment_name=experiment_name,
        average=user_defined_mfa_settings.average_data,
        dataset_id=dataset_id, **user_defined_mfa_settings.data_content_parameter_dict)
    user_defined_model = data_model_object.user_defined_model
    mfa_model_obj = data_model_object.mfa_model_obj
    solver_parameter_dict = user_defined_mfa_settings.construct_solver_parameter_dict(
        test_mode=test_mode, docker_mode=docker_mode, suffix=suffix, parallel_num=parallel_num)
    slsqp_mfa_config = user_defined_mfa_settings.construct_mfa_config({
        'user_defined_specific_flux_range_dict': user_defined_mfa_settings.user_defined_specific_flux_range_dict,
        'user_defined_constant_flux_dict': user_defined_mfa_settings.user_defined_constant_flux_dict,
        'model_tissue_set': user_defined_model.model_tissue_set,
        'flux_name_index_dict': mfa_model_obj.flux_name_index_dict})

    # 构建 parameter_label_content_dict
    parameter_label_content_dict = {}
    for result_label, mfa_data_obj in mfa_data_dict.items():
        parameter_label_content_dict[result_label] = (
            (None, result_label, None),
            (mfa_model_obj, mfa_data_obj, slsqp_mfa_config),
            result_information_dict[result_label],
            {}
        )
    return (
        parameter_label_content_dict, solver_parameter_dict, mfa_data_dict,
        experimental_result_information_dict_analysis,
        experimental_tissue_mid_name_list_dict_constructor, experimental_figure_config_dict,
        result_information_dict)


def model_data_config_initialization(
        experiment_name, running_mode, test_mode, docker_mode, dataset_id, suffix, parallel_num):
    from ..model_data_combination.common_data_model_loader import common_data_model_function_loader
    data_model_object = common_data_model_function_loader(experiment_name)

    user_defined_mfa_settings = data_model_object.user_defined_mfa_settings
    user_defined_mfa_settings.initialize_running_mode(running_mode, test_mode)

    # 尝试获取多模型映射（新API）；若不存在则回退到旧的单模型方式
    try:
        (
            result_label_model_data_config_dict, experimental_result_information_dict_analysis,
            experimental_tissue_mid_name_list_dict_constructor,
            result_information_dict, experimental_figure_config_dict
        ) = (
            data_model_object.return_result_label_model_data_config_dict(
                running_mode=running_mode, experiment_name=experiment_name,
                dataset_id=dataset_id, suffix=suffix))
    except AttributeError as e:
        (
            parameter_label_content_dict, solver_parameter_dict,
            mfa_data_dict_for_experimental_plotting, experimental_result_information_dict_analysis,
            experimental_tissue_mid_name_list_dict_constructor, experimental_figure_config_dict,
            result_information_dict
        ) = old_model_data_config_generator(
            data_model_object, user_defined_mfa_settings=user_defined_mfa_settings,
            running_mode=running_mode, experiment_name=experiment_name,
            dataset_id=dataset_id, test_mode=test_mode, docker_mode=docker_mode, suffix=suffix,
            parallel_num=parallel_num)
    else:
        # 统一流程：按 result_label 生成各自的 config
        parameter_label_content_dict = {}
        mfa_data_dict_for_experimental_plotting = {}
        solver_parameter_dict = user_defined_mfa_settings.construct_solver_parameter_dict(
            test_mode=test_mode, docker_mode=docker_mode, suffix=suffix, parallel_num=parallel_num)
        for result_label, (
                model_label, mfa_model_obj,
                data_label, mfa_data_obj,
                config_label, slsqp_mfa_config,
        ) in result_label_model_data_config_dict.items():
            parameter_label_content_dict[result_label] = (
                (model_label, data_label, config_label),
                (mfa_model_obj, mfa_data_obj, slsqp_mfa_config),
                result_information_dict[result_label],
                {}
            )
            mfa_data_dict_for_experimental_plotting[result_label] = mfa_data_obj

    if user_defined_mfa_settings.tf2_solver:
        result_class = OutOfOrderFinalResult
    else:
        result_class = FinalResult

    result_process_information_dict_analysis = data_model_object.result_process_information_dict_analysis
    try:
        predicted_tissue_mid_name_list_dict_constructor = data_model_object.tissue_mid_name_list_dict_constructor
    except AttributeError:
        predicted_tissue_mid_name_list_dict_constructor = None
    try:
        tissue_flux_name_list_dict_constructor = data_model_object.tissue_flux_name_list_dict_constructor
    except AttributeError:
        tissue_flux_name_list_dict_constructor = None
    try:
        result_process_figure_config_dict = data_model_object.result_process_figure_config_dict_generator()
    except AttributeError:
        result_process_figure_config_dict = {}

    experimental_data_display_parameter_dict = {
        'result_information_dict_analysis': experimental_result_information_dict_analysis,
        'tissue_mid_name_list_dict_constructor': experimental_tissue_mid_name_list_dict_constructor,
        'experimental_figure_config_dict': experimental_figure_config_dict,
    }
    final_process_parameter_dict = {
        'result_process_information_dict_analysis': result_process_information_dict_analysis,
        'tissue_mid_name_list_dict_constructor': predicted_tissue_mid_name_list_dict_constructor,
        'tissue_flux_name_list_dict_constructor': tissue_flux_name_list_dict_constructor,
        **user_defined_mfa_settings.final_process_parameter_dict,
        **result_process_figure_config_dict,
    }
    return (
        parameter_label_content_dict, mfa_data_dict_for_experimental_plotting,
        user_defined_mfa_settings, result_information_dict,
        solver_parameter_dict, experimental_data_display_parameter_dict, final_process_parameter_dict, result_class,
    )


def multiple_tissue_model_analysis_main(
        running_mode, experiment_name=None, dataset_id=None, suffix=None, test_mode=False, docker_mode=None,
        parallel_num=None):
    if experiment_name is None:
        experiment_name = default_data_model_name
    (
        parameter_label_content_dict, mfa_data_dict_for_experimental_plotting,
        user_defined_mfa_settings, result_information_dict,
        solver_parameter_dict, experimental_data_display_parameter_dict,
        final_process_parameter_dict, result_class
    ) = model_data_config_initialization(
        experiment_name, running_mode, test_mode, docker_mode, dataset_id, suffix, parallel_num)

    final_result_obj = result_class(
        Direct.output_direct, Direct.common_data_direct, experiment_name, suffix=suffix)

    if running_mode == MFARunningMode.flux_analysis:
        exit_code = common_solver(
            parameter_label_content_dict, final_result_obj, **solver_parameter_dict)
        sys.exit(exit_code)
    else:
        if running_mode == MFARunningMode.raw_experimental_data_plotting:
            experimental_mid_and_raw_data_plotting(
                mfa_data_dict_for_experimental_plotting, result_information_dict, final_result_obj, **experimental_data_display_parameter_dict)
        else:
            solver_dict, same_model_dict, same_data_dict = solver_dict_constructor(parameter_label_content_dict)
            if running_mode == MFARunningMode.solver_output:
                solver_output(
                    solver_dict, result_information_dict, final_result_obj, same_model_dict, same_data_dict)
            elif running_mode == MFARunningMode.result_process:
                final_result_obj.final_process(
                    solver_dict=solver_dict, **final_process_parameter_dict)
            else:
                raise ValueError()

