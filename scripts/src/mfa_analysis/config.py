from ..core.common.classes import OptionDict, MFAConfig, TransformDict
from ..core.common.functions import (
    tissue_specific_name_constructor, check_if_mix_flux, mid_name_process, tissue_name_breakdown)
from ..core.common.config import ParamName, ModelKeyword, CoreConstants
from ..common_parallel_solver.config import Keywords as ParallelSolverKeywords
from ..common_parallel_solver.result_class import FinalResult as GeneralFinalResult
from ..io_functions.config import check_and_mkdir_of_direct
from ..io_functions.result_processing_functions import (
    multiple_repeat_result_processing, loss_data_distribution_plotting,
    time_distribution_plotting, experimental_mid_prediction, common_flux_comparison_func, figure_plotting,
    experimental_data_plotting_func_template
)
from ..io_functions.result_output_functions import output_raw_flux_data, output_predicted_mid_data

from ..common.config import (
    Direct as CommonDirect, Keywords as GeneralKeywords, Color, FigureData, FigureDataKeywords,
)
from ..common.packages import mp, np, os
from ..common.functions import (
    update_parameter_object, reversible_flux_pair_dict_generator,
    tissue_specific_reversible_flux_title_constructor_generator,)
from scripts.model.model_loader import ModelList, model_loader
from scripts.data.data_loader import DataSource, common_data_loader

from ..inventory import MFARunningMode, ExperimentName, default_data_model_name
from ..model_data_combination.common_data_model_loader import (
    common_data_model_function_loader)


class Direct(object):
    name = 'mfa_analysis'
    output_direct = f'{CommonDirect.output_direct}/{name}'
    common_data_direct = f'{CommonDirect.common_submitted_raw_data_direct}/{name}'
    experimental_result_display = 'experimental_result_display'


class Keywords(GeneralKeywords):
    pass


# test_mode = True
# running_mode = MFARunningMode.flux_analysis
# running_mode = MFARunningMode.result_process
running_mode = MFARunningMode.raw_experimental_data_plotting


# def running_settings(
#         user_defined_mfa_settings, experiment_name, running_mode, test_mode, docker_mode, suffix=None, parallel_num=None,
#         average_data=False, benchmark=False):
#     raw_slsqp_mfa_config = user_defined_mfa_settings.slsqp_mfa_config
#
#     tf2_benchmark_solution_num = 1000
#     # tf2_benchmark_solution_num = 200
#     slsqp_benchmark_solution_num = 100
#     batch_size_32g = 90
#     batch_size_64g = 150
#     batch_size_96g = 240
#     nnls_parallel_num_32g_64g = 30
#     nnls_parallel_num_96g = 40
#     if running_mode == MFARunningMode.flux_analysis_tf2:
#         mp.set_start_method('spawn')
#         # if os.name == 'nt':
#         #     mp.set_start_method('spawn')
#         # else:
#         #     mp.set_start_method('forkserver')
#         parameter_tuple_for_each_suffix_dict = {
#             None: (batch_size_32g, 4000, nnls_parallel_num_32g_64g),
#             '0': (batch_size_32g, 4000, nnls_parallel_num_32g_64g),
#             '1': (batch_size_64g, 5000, nnls_parallel_num_32g_64g),
#             '2': (batch_size_96g, 4000, nnls_parallel_num_96g),
#         }
#         result_label_obj_value_dict = {
#             '1599': 0.7, '1571': 1.3, '1593': 1.1,
#         }
#         if test_mode:
#             batch_size = 39
#             nnls_batch_size = 3
#             nnls_parallel_num = batch_size
#             emu_value_parallel_num = nnls_parallel_num
#             each_portion_flux_num = 2000
#             max_optimization_each_generation = batch_size
#             emu_gradient_hessian_parallel_num = 3
#         else:
#             emu_gradient_hessian_parallel_num = 3
#             (batch_size, each_portion_flux_num, nnls_parallel_num) = parameter_tuple_for_each_suffix_dict[suffix]
#             # nnls_parallel_num = 30
#             # emu_value_parallel_num = 30
#             emu_value_parallel_num = nnls_parallel_num
#             max_optimization_each_generation = int(0.35 * batch_size)
#         report_interval = 0
#         if benchmark:
#             if test_mode:
#                 each_case_optimization_num = 2 * batch_size
#             else:
#                 each_case_optimization_num = tf2_benchmark_solution_num
#             slsqp_mfa_config = {
#                 result_label: user_defined_mfa_settings.construct_tf2_slsqp_mfa_config(
#                     raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num,
#                     emu_gradient_hessian_parallel_num, target_obj_value)
#                 for result_label, target_obj_value in result_label_obj_value_dict.items()}
#         else:
#             each_case_optimization_num = 800
#             slsqp_mfa_config = user_defined_mfa_settings.construct_tf2_slsqp_mfa_config(
#                 raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num,
#                 emu_gradient_hessian_parallel_num)
#         parallel_parameter_dict = {
#             Keywords.max_optimization_each_generation: max_optimization_each_generation,
#             Keywords.batch_solving: batch_size,
#         }
#         print(f'Current experiment: Batch size: {batch_size}, Hessian portion size: {each_portion_flux_num}, '
#               f'EMU gradient/hessian parallel num: {emu_gradient_hessian_parallel_num}, '
#               f'NNLS parallel num: {nnls_parallel_num}, EMU value parallel num: {emu_value_parallel_num}')
#     else:
#         mp.set_start_method('spawn')
#         report_interval = 1
#         parallel_parameter_dict = {
#             Keywords.max_optimization_each_generation: 10,
#             Keywords.each_process_optimization_num: 5,
#             Keywords.processes_num: 6,
#             Keywords.thread_num_constraint: 4,
#             # Keywords.parallel_test: True,
#         }
#         if average_data:
#             optimization_num_for_each_suffix_dict = {
#                 None: 900,
#                 '0': 700,
#                 '1': 600,
#                 '2': 200,
#             }
#         else:
#             optimization_num_for_each_suffix_dict = {
#                 None: 280,
#                 '0': 100,
#                 '1': 20,
#                 '2': 20,
#             }
#         if test_mode:
#             each_case_optimization_num = 2
#             parallel_parameter_dict = None
#         elif docker_mode:
#             each_case_optimization_num = 20
#             parallel_parameter_dict = None
#         else:
#             if benchmark:
#                 each_case_optimization_num = slsqp_benchmark_solution_num
#             else:
#                 each_case_optimization_num = optimization_num_for_each_suffix_dict[suffix]
#             if parallel_num is not None:
#                 assert isinstance(parallel_num, int)
#                 parallel_parameter_dict[Keywords.processes_num] = parallel_num
#         slsqp_mfa_config = raw_slsqp_mfa_config
#     if 'tf2' in experiment_name.value:
#         result_parallel_num = batch_size_32g
#     else:
#         result_parallel_num = 7
#     if benchmark:
#         average_optimized_results = 3
#         loss_percentile = 9 / 200
#         result_label_suffix_tuple = tuple()
#         user_defined_mfa_settings.update_final_process_parameter_dict(
#             benchmark=benchmark, parallel_num=result_parallel_num, average_optimized_results=average_optimized_results,
#             loss_percentile=loss_percentile, result_label_suffix_tuple=result_label_suffix_tuple)
#     else:
#         user_defined_mfa_settings.update_final_process_parameter_dict(
#             benchmark=benchmark, parallel_num=1)
#     return (
#         report_interval, each_case_optimization_num, parallel_parameter_dict, slsqp_mfa_config)


