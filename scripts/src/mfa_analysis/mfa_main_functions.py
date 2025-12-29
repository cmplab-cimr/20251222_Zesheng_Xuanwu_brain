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


def model_data_config_initialization(
        experiment_name, running_mode, test_mode, docker_mode, dataset_id, suffix, parallel_num):
    from ..model_data_combination.common_data_model_loader import common_data_model_function_loader
    data_model_object = common_data_model_function_loader(experiment_name)

    data_name = data_model_object.data_name
    multi_tissue_user_defined_model = data_model_object.user_defined_model
    multi_organ_mfa_model = data_model_object.mfa_model_obj
    user_defined_mfa_settings = data_model_object.user_defined_mfa_settings
    user_defined_mfa_settings.initialize_running_mode(running_mode, test_mode)
    user_defined_mfa_settings.specific_flux_range_constant_flux_dict_modifier(
        model_tissue_set=multi_tissue_user_defined_model.model_tissue_set,
        flux_name_index_dict=multi_organ_mfa_model.flux_name_index_dict)
    average_data = user_defined_mfa_settings.average_data

    if user_defined_mfa_settings.tf2_solver:
        result_class = OutOfOrderFinalResult
    else:
        result_class = FinalResult

    (
        mfa_data_dict, experimental_result_information_dict_analysis, experimental_tissue_mid_name_list_dict_constructor,
        result_information_dict, experimental_figure_config_dict,
    ) = data_model_object.return_data_content(
        running_mode=running_mode, experiment_name=experiment_name, average=average_data,
        dataset_id=dataset_id, **user_defined_mfa_settings.data_content_parameter_dict)

    keyword = data_model_object.keyword
    result_process_information_dict_analysis = data_model_object.result_process_information_dict_analysis
    predicted_tissue_mid_name_list_dict_constructor = data_model_object.tissue_mid_name_list_dict_constructor
    tissue_flux_name_list_dict_constructor = data_model_object.tissue_flux_name_list_dict_constructor
    result_process_figure_config_dict = data_model_object.result_process_figure_config_dict_generator()

    (slsqp_mfa_config, solver_parameter_dict) = user_defined_mfa_settings.construct_running_settings(
        test_mode, docker_mode, suffix, parallel_num,
    )
    experimental_data_display_parameter_dict = {
        'result_information_dict_analysis': experimental_result_information_dict_analysis,
        'tissue_mid_name_list_dict_constructor': experimental_tissue_mid_name_list_dict_constructor,
        'experimental_figure_config_dict': experimental_figure_config_dict,
        'average_data': average_data,
    }
    final_process_parameter_dict = {
        'result_process_information_dict_analysis': result_process_information_dict_analysis,
        'tissue_mid_name_list_dict_constructor': predicted_tissue_mid_name_list_dict_constructor,
        'tissue_flux_name_list_dict_constructor': tissue_flux_name_list_dict_constructor,
        **user_defined_mfa_settings.final_process_parameter_dict,
        **result_process_figure_config_dict,
    }
    return (
        multi_organ_mfa_model, mfa_data_dict, user_defined_mfa_settings, result_information_dict, slsqp_mfa_config,
        solver_parameter_dict, experimental_data_display_parameter_dict, final_process_parameter_dict, result_class,
    )


def multiple_tissue_model_analysis_main(
        running_mode, experiment_name=None, dataset_id=None, suffix=None, test_mode=False, docker_mode=None,
        parallel_num=None):
    if experiment_name is None:
        experiment_name = default_data_model_name
    (
        multi_organ_mfa_model, mfa_data_dict, user_defined_mfa_settings, result_information_dict,
        slsqp_mfa_config, solver_parameter_dict, experimental_data_display_parameter_dict,
        final_process_parameter_dict, result_class
    ) = model_data_config_initialization(
        experiment_name, running_mode, test_mode, docker_mode, dataset_id, suffix, parallel_num)

    final_result_obj = result_class(
        Direct.output_direct, Direct.common_data_direct, experiment_name, suffix=suffix)

    parameter_label_content_dict = {}
    for result_label, mfa_data_obj in mfa_data_dict.items():
        data_label = result_label
        result_information = result_information_dict[result_label]
        other_information_dict = {}
        if isinstance(slsqp_mfa_config, dict):
            new_mfa_config = slsqp_mfa_config[result_label].copy()
        else:
            new_mfa_config = slsqp_mfa_config.copy()
        parameter_label_content_dict[result_label] = (
            (None, data_label, None), (multi_organ_mfa_model, mfa_data_obj, new_mfa_config),
            result_information, other_information_dict
        )

    if running_mode == MFARunningMode.flux_analysis:
        exit_code = common_solver(
            parameter_label_content_dict, final_result_obj, **solver_parameter_dict)
        sys.exit(exit_code)
    else:
        if running_mode == MFARunningMode.raw_experimental_data_plotting:
            experimental_mid_and_raw_data_plotting(
                mfa_data_dict, result_information_dict, final_result_obj, **experimental_data_display_parameter_dict)
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

