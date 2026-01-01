from ..core.common.classes import OptionDict, MFAConfig, TransformDict
from ..core.common.functions import tissue_specific_name_constructor, check_if_mix_flux
from ..core.common.config import ParamName, ModelKeyword, CoreConstants
from ..common_parallel_solver.packages import mp
from ..common_parallel_solver.config import Keywords as ParallelSolverKeywords

from ..inventory import MFARunningMode


class MFASettings(object):
    tf2_solver = False

    loss_percentile = 0.005
    report_interval = 50
    thread_num_constraint = None
    update_data = False
    solver_type = ParamName.slsqp_numba_python_solver
    load_previous_results = True
    average_optimized_results = 2
    tf2_batch_size = 10
    extra_result_suffix_tuple = (None, 0, 1, 2, 3)
    average_data = False
    benchmark = False

    optimization_num_for_each_suffix_dict = {
        None: 0
    }

    data_content_parameter_dict = {}

    user_defined_specific_flux_range_dict = {
        'Salvage_c': (None, (1, 10)),
        'G6PDH2R_PGL_GND_c': (None, (1, 200)),
    }

    user_defined_constant_flux_dict = {
        'GLC_infuse': (ModelKeyword.serum, 100)
    }

    common_flux_range = (1, 3000)
    common_mix_ratio_range = (0.05, 0.95)
    mix_ratio_multiplier = 100

    # Parameters in slsqp_solver_config_dict
    loss_type = ParamName.mean_squared_loss
    debug = True
    slsqp_max_iter = 2000

    final_process_parameter_dict = {}

    result_label_obj_value_dict = None

    initial_slsqp_multiple = 2
    initial_batch_multiple = 10
    initial_processes_num = 20

    def __init__(self):
        self.slsqp_solver_config_dict = OptionDict({
            ParamName.loss_type: self.loss_type,
            ParamName.debug: self.debug,
            ParamName.slsqp_max_iter: self.slsqp_max_iter,
        })
        self.default_slsqp_mfa_config = MFAConfig(
            common_flux_range=self.common_flux_range,
            specific_flux_range_dict=None, dynamic_constant_flux_list=[],
            preset_constant_flux_value_dict=None,
            common_mix_ratio_range=self.common_mix_ratio_range, mix_ratio_multiplier=self.mix_ratio_multiplier,
            solver_type=self.solver_type,
            solver_config_dict=self.slsqp_solver_config_dict)

    def initialize_running_mode(self, running_mode, test_mode, **kwargs):
        if running_mode != MFARunningMode.flux_analysis or test_mode:
            self.tf2_solver = False

    @staticmethod
    def construct_tf2_slsqp_mfa_config(
            raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num,
            emu_gradient_hessian_parallel_num, target_loss_value=None):
        if target_loss_value is None:
            target_loss_value = 99999
        tf2_slsqp_solver_config_dict = OptionDict({
            ParamName.batch_size: batch_size,
            ParamName.each_portion_flux_num: each_portion_flux_num,
            ParamName.refresh_hessian_threshold: 3,
            ParamName.less_update_ratio_for_hessian_refresh: 0.05,
            ParamName.sqp_less_update_threshold_ratio: 0.5,
            ParamName.hessian_reset_threshold_time: 10,
            # ParamName.nnls_batch_size: nnls_batch_size,
            ParamName.nnls_parallel_num: nnls_parallel_num,
            ParamName.hessian_refresh_block_time: 6,
            ParamName.hessian_refresh_max_interval: 15,
            ParamName.recent_mid_monitor_num: 10,
            ParamName.recent_mid_average_update_threshold: 5e-5,
            ParamName.recent_average_flux_update_threshold: 2,
            ParamName.final_loss_threshold: target_loss_value,
            # ParamName.emu_prediction_parallel_num: emu_prediction_parallel_num,
            ParamName.emu_gradient_hessian_parallel_num: emu_gradient_hessian_parallel_num,
            ParamName.emu_value_parallel_num: emu_value_parallel_num,
        })
        tf2_slsqp_mfa_config = raw_slsqp_mfa_config.copy()
        tf2_slsqp_mfa_config.solver_type = ParamName.tf2_slsqp_solver
        tf2_slsqp_mfa_config.solver_config_dict = tf2_slsqp_solver_config_dict
        return tf2_slsqp_mfa_config

    def specific_flux_range_constant_flux_dict_modifier(self, model_tissue_set, flux_name_index_dict,):
        def add_new_component(_flux_dict, _flux_name, _flux_tissue, _flux_range):
            _new_flux_name = tissue_specific_name_constructor(_flux_name, tissue_name_list=_flux_tissue)
            if _new_flux_name in flux_name_index_dict:
                _flux_dict[_new_flux_name] = _flux_range

        def construct_flux_dict(user_defined_flux_dict):
            target_flux_dict = {}
            for flux_name, (flux_tissue, flux_range) in user_defined_flux_dict.items():
                if flux_tissue is not None:
                    if isinstance(flux_tissue, str):
                        add_new_component(target_flux_dict, flux_name, flux_tissue, flux_range)
                    elif isinstance(flux_tissue, (tuple, list)):
                        for flux_tissue_item in flux_tissue:
                            add_new_component(target_flux_dict, flux_name, flux_tissue_item, flux_range)
                    else:
                        raise ValueError('Type not recognized!')
                else:
                    for model_tissue in model_tissue_set:
                        add_new_component(target_flux_dict, flux_name, model_tissue, flux_range)
            return target_flux_dict

        specific_flux_range_dict = construct_flux_dict(self.user_defined_specific_flux_range_dict)
        constant_flux_dict = construct_flux_dict(self.user_defined_constant_flux_dict)
        raw_mfa_config = self.default_slsqp_mfa_config
        if isinstance(raw_mfa_config, dict):
            for each_tissue_raw_mfa_config in raw_mfa_config.values():
                each_tissue_raw_mfa_config.specific_flux_range_dict = dict(specific_flux_range_dict)
                each_tissue_raw_mfa_config.preset_constant_flux_value_dict = dict(constant_flux_dict)
        else:
            raw_mfa_config.specific_flux_range_dict = specific_flux_range_dict
            raw_mfa_config.preset_constant_flux_value_dict = constant_flux_dict

    def _generate_solver_parameter_dict(
            self, each_case_target_optimization_num, test_mode=False, report_interval=None,
            parallel_parameter_dict=None, load_previous_results=False, docker_mode=False, **kwargs):
        solver_parameter_dict = {
            'each_case_target_optimization_num': each_case_target_optimization_num,
            'test_mode': test_mode,
            'parallel_parameter_dict': parallel_parameter_dict,
            'load_results': load_previous_results,
            'docker_mode': docker_mode,
            'report_interval': report_interval,
        }
        solver_parameter_dict.update(kwargs)
        return solver_parameter_dict

    def update_final_process_parameter_dict(
            self, result_parallel_num=None, average_optimized_results=None,
            loss_percentile=None, result_label_suffix_tuple=None, **kwargs):
        if average_optimized_results is None:
            average_optimized_results = self.average_optimized_results
        if loss_percentile is None:
            loss_percentile = self.loss_percentile
        if result_label_suffix_tuple is None:
            result_label_suffix_tuple = self.extra_result_suffix_tuple
        self.final_process_parameter_dict = {
            'parallel_num': result_parallel_num,
            'average_optimized_results': average_optimized_results,
            'result_label_suffix_tuple': result_label_suffix_tuple,
            'loss_percentile': loss_percentile,
            **kwargs,
        }

    def tf2_running_parameters(self):
        batch_size_32g = 90
        batch_size_64g = 150
        batch_size_96g = 240
        nnls_parallel_num_32g_64g = 30
        nnls_parallel_num_96g = 40
        test_batch_size = 39
        parameter_tuple_for_each_suffix_dict = {
            'test': (test_batch_size, 2000, test_batch_size, test_batch_size),
            None: (batch_size_32g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g),
            '0': (batch_size_32g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g),
            '1': (batch_size_64g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g),
            '2': (batch_size_96g, 4000, nnls_parallel_num_96g, nnls_parallel_num_96g),
        }
        return batch_size_32g, parameter_tuple_for_each_suffix_dict

    def _construct_tf2_running_settings(
            self, suffix, test_mode, docker_mode, parallel_num):
        raw_slsqp_mfa_config = self.default_slsqp_mfa_config

        mp.set_start_method('spawn')
        # if os.name == 'nt':
        #     mp.set_start_method('spawn')
        # else:
        #     mp.set_start_method('forkserver')
        result_parallel_num, parameter_tuple_for_each_suffix_dict = self.tf2_running_parameters()
        if test_mode:
            (
                batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num
            ) = parameter_tuple_for_each_suffix_dict['test']
            each_case_target_optimization_num = 2 * batch_size
        else:
            (
                batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num
            ) = parameter_tuple_for_each_suffix_dict[suffix]
            each_case_target_optimization_num = self.optimization_num_for_each_suffix_dict[suffix]
        max_optimization_each_generation = self.initial_batch_multiple * batch_size
        initial_processes_num = self.initial_processes_num
        emu_gradient_hessian_parallel_num = 3

        if self.result_label_obj_value_dict is not None:
            slsqp_mfa_config = {
                result_label: self.construct_tf2_slsqp_mfa_config(
                    raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num,
                    emu_value_parallel_num,
                    emu_gradient_hessian_parallel_num, target_obj_value)
                for result_label, target_obj_value in self.result_label_obj_value_dict.items()}
        else:
            slsqp_mfa_config = self.construct_tf2_slsqp_mfa_config(
                raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num, emu_value_parallel_num,
                emu_gradient_hessian_parallel_num)
        parallel_parameter_dict = {
            ParallelSolverKeywords.max_optimization_each_generation: max_optimization_each_generation,
            ParallelSolverKeywords.batch_solving: batch_size,
            ParallelSolverKeywords.processes_num: initial_processes_num,
        }
        self.update_final_process_parameter_dict(result_parallel_num=result_parallel_num)
        solver_parameter_dict = self._generate_solver_parameter_dict(
            each_case_target_optimization_num=each_case_target_optimization_num,
            test_mode=test_mode, parallel_parameter_dict=parallel_parameter_dict,
            load_previous_results=self.load_previous_results, docker_mode=docker_mode)
        print(f'Current experiment: suffix: {suffix}, Batch size: {batch_size}, Hessian portion size: {each_portion_flux_num}, '
              f'EMU gradient/hessian parallel num: {emu_gradient_hessian_parallel_num}, '
              f'NNLS parallel num: {nnls_parallel_num}, EMU value parallel num: {emu_value_parallel_num}')
        return slsqp_mfa_config, solver_parameter_dict

    def ordinary_slsqp_running_parameter(self):
        each_process_optimization_num = 20
        report_interval = int(each_process_optimization_num / 5)
        return each_process_optimization_num, report_interval

    def _construct_ordinary_slsqp_running_settings(
            self, suffix, test_mode, docker_mode, parallel_num=None):
        mp.set_start_method('spawn')
        each_process_optimization_num, report_interval = self.ordinary_slsqp_running_parameter()

        if test_mode:
            each_case_target_optimization_num = 2
            parallel_parameter_dict = None
        elif docker_mode:
            each_case_target_optimization_num = int(self.optimization_num_for_each_suffix_dict[suffix] / parallel_num)
            parallel_parameter_dict = None
        else:
            each_case_target_optimization_num = self.optimization_num_for_each_suffix_dict[suffix]
            default_parallel_num = 6
            if parallel_num is None:
                parallel_num = default_parallel_num
            else:
                assert isinstance(parallel_num, int)
            max_optimization_each_generation = min(
                int(each_process_optimization_num * parallel_num * self.initial_slsqp_multiple),
                each_case_target_optimization_num)
            parallel_parameter_dict = {
                ParallelSolverKeywords.max_optimization_each_generation: max_optimization_each_generation,
                ParallelSolverKeywords.each_process_optimization_num: each_process_optimization_num,
                ParallelSolverKeywords.processes_num: parallel_num,
                ParallelSolverKeywords.thread_num_constraint: 4,
            }
            print(f'Current experiment: suffix: {suffix}, max initial value: {max_optimization_each_generation}, '
                  f'each process max num: {each_process_optimization_num}, parallel num: {parallel_num}')
        slsqp_mfa_config = self.default_slsqp_mfa_config
        result_parallel_num = 7
        self.update_final_process_parameter_dict(
            result_parallel_num=result_parallel_num, average_optimized_results=None,
            loss_percentile=None, result_label_suffix_tuple=None)
        solver_parameter_dict = self._generate_solver_parameter_dict(
            each_case_target_optimization_num=each_case_target_optimization_num,
            test_mode=test_mode, report_interval=report_interval,
            parallel_parameter_dict=parallel_parameter_dict,
            load_previous_results=self.load_previous_results, docker_mode=docker_mode)
        return slsqp_mfa_config, solver_parameter_dict

    def construct_running_settings(
            self, test_mode, docker_mode, suffix, parallel_num):
        if self.tf2_solver:
            slsqp_mfa_config, solver_parameter_dict = self._construct_tf2_running_settings(
                suffix, test_mode, docker_mode, parallel_num)
        else:
            slsqp_mfa_config, solver_parameter_dict = self._construct_ordinary_slsqp_running_settings(
                suffix, test_mode, docker_mode, parallel_num)

        return slsqp_mfa_config, solver_parameter_dict

# benchmark parameter:

# average_optimized_results = 3
# loss_percentile = 9 / 200
# result_label_suffix_tuple = tuple()
# slsqp_benchmark_solution_num = 100
# each_case_target_optimization_num = slsqp_benchmark_solution_num
# if test_mode:
#     each_case_target_optimization_num = 2 * batch_size
# else:
#     each_case_target_optimization_num = tf2_benchmark_solution_num
# slsqp_mfa_config = {
#     result_label: self.construct_tf2_slsqp_mfa_config(
#         raw_slsqp_mfa_config, batch_size, each_portion_flux_num, nnls_parallel_num,
#         emu_value_parallel_num,
#         emu_gradient_hessian_parallel_num, target_obj_value)
#     for result_label, target_obj_value in result_label_obj_value_dict.items()}
# tf2_benchmark_solution_num = 1000
# tf2_benchmark_solution_num = 200
# result_label_obj_value_dict = {
#     '1599': 0.7, '1571': 1.3, '1593': 1.1,
# }

