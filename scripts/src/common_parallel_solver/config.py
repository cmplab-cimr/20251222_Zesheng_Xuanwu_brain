from .packages import np, gzip, pickle, pathlib

from ..core.common.functions import ProgressBarExtraProcess
from ..core.solver.solver_construction_functions.solver_constructor import specific_solver_constructor, \
    base_solver_constructor, common_solver_constructor
from ..core.sampler.np_sampler.sampler_class import OptGpSampler

random_seed = np.random.default_rng(4536251)


class Keywords(object):
    thread_num_constraint = 'thread_num_constraint'
    each_process_optimization_num = 'each_process_optimization_num'
    max_optimization_each_generation = 'max_optimization_each_generation'
    specific_target_optimization_num = 'specific_target_optimization_num'
    predefined_initial_solution_matrix = 'predefined_initial_solution_matrix'
    parallel_test = 'parallel_test'
    processes_num = 'processes_num'
    unoptimized = 'unoptimized'
    batch_solving = 'batch_solving'
    maximal_save_point = 'maximal_save_point'

    flux_raw_data = 'flux_raw_data'
    mid_raw_data = 'mid_raw_data'
    solver_descriptions = 'solver_descriptions'
    model_metabolites_reactions_standard_name = 'model_metabolites_reactions_standard_name'


class Direct(object):
    raw_flux_analysis = 'raw_flux_analysis'
    flux_comparison = 'flux_comparison'
    experimental_result_display = 'experimental_result_display'
    predicted_experimental_mid_comparison_direct = 'predicted_experimental_mid_comparison'
    metabolic_network_visualization_direct = 'metabolic_network_visualization'
    flux_result_xlsx_filename = f'{Keywords.flux_raw_data}.xlsx'
    mid_result_xlsx_filename = f'{Keywords.mid_raw_data}.xlsx'
    solver_description_xlsx_filename = f'{Keywords.solver_descriptions}.xlsx'

    solution_array = 'solution_array'
    time_array = 'time_array'
    loss_array = 'loss_array'
    result_information = 'result_information'
    predicted_dict = 'predicted_dict'
    flux_name_index_dict = 'flux_name_index_dict'
    experimental_data = 'experimental_data'
    solution_id_array = 'solution_id_array'


def parameter_extract(parameter_dict, parameter_name, default_value=None):
    if parameter_name in parameter_dict:
        return parameter_dict[parameter_name]
    else:
        return default_value


def pickle_save(obj, file_path):
    with gzip.open(file_path, 'wb') as f_out:
        pickle.dump(obj, f_out)


def pickle_load(file_path):
    with gzip.open(file_path, 'rb') as f_in:
        obj = pickle.load(f_in)
    return obj


def check_and_mkdir_of_direct(direct_str, file_path=False):
    direct_obj = pathlib.Path(direct_str)
    if file_path:
        direct_obj = direct_obj.parent
    dir_stack = []
    while not direct_obj.exists():
        dir_stack.append(direct_obj)
        direct_obj = direct_obj.parent
    while len(dir_stack) != 0:
        missed_direct = dir_stack.pop()
        missed_direct.mkdir()


def npz_load(raw_path, *args, allow_pickle=False):
    result_list = []
    suffix_path = '{}.npz'.format(raw_path)
    with np.load(suffix_path, allow_pickle=allow_pickle) as data:
        for label in args:
            result_list.append(data[label])
    if len(result_list) == 1:
        return result_list[0]
    else:
        return result_list


def npz_save(path, **kwargs):
    np.savez_compressed(path, **kwargs)


def split_total_num_to_process(total_num, processes_num):
    common_size = total_num // processes_num
    rest_size = total_num % processes_num
    each_process_num_list = (
            [common_size + 1] * rest_size +
            [common_size] * (processes_num - rest_size))
    return each_process_num_list
