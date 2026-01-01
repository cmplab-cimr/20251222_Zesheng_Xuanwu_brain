from ..common.config import CoreConstants
from ..common.packages import np, optimize, mp
from ..common.functions import split_total_num_to_process, ProgressBarExtraProcess

default_random_min_factor = CoreConstants.lp_random_min_factor
default_random_max_factor = CoreConstants.lp_random_max_factor
default_maximal_failed_time = CoreConstants.maximal_failed_time
eps = CoreConstants.eps_for_computation


class BaseSampler(object):
    def __init__(
            self, variable_num, eq_matrix=None, eq_right_side_vector=None,
            smaller_eq_matrix=None, smaller_eq_right_side_vector=None,
            min_value_vector=None, max_value_vector=None, dynamic_constraint_index_dict=None,
            numeric_eps=eps, verbose=False):
        if len(min_value_vector) != variable_num:
            raise ValueError("Size of min_value_vector should equal to variable number")
        if len(max_value_vector) != variable_num:
            raise ValueError("Size of max_value_vector should equal to variable number")
        if smaller_eq_matrix is not None and smaller_eq_right_side_vector is not None:
            if smaller_eq_matrix.shape[1] != variable_num:
                raise ValueError("Column of smaller eq matrix should equal to variable number")
            if smaller_eq_matrix.shape[0] != smaller_eq_right_side_vector.shape[0]:
                raise ValueError("Row of ineq matrix should equal to right side")
        elif smaller_eq_matrix is not None or smaller_eq_right_side_vector is not None:
            raise ValueError("")
        if eq_matrix is not None and eq_right_side_vector is not None:
            if eq_matrix.shape[1] != variable_num:
                raise ValueError("Upper bound should equal to variable number")
            if eq_matrix.shape[0] != eq_right_side_vector.shape[0]:
                raise ValueError("Row of ineq matrix should equal to right side")
        elif eq_matrix is not None or eq_right_side_vector is not None:
            raise ValueError("")
        self.variable_num = variable_num
        self.eq_matrix = eq_matrix
        self.eq_right_side_vector = eq_right_side_vector
        self.smaller_eq_matrix = smaller_eq_matrix
        self.smaller_eq_right_side_vector = smaller_eq_right_side_vector
        self.min_value_vector = min_value_vector
        self.max_value_vector = max_value_vector
        if numeric_eps is None:
            numeric_eps = CoreConstants.eps_for_sampling
        self.numeric_eps = numeric_eps
        self.verbose = verbose
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        self.dynamic_constraint_index_dict = dynamic_constraint_index_dict
        self._set_dynamic_constraint = len(dynamic_constraint_index_dict) == 0

    def refresh(self):
        """Abstract refresh function.

        Should be overwritten by child classes.

        """
        pass

    def set_right_side_vector(self, dynamic_flux_value_dict):
        self._set_dynamic_constraint = True

    def sample(self, n=1):
        if not self._set_dynamic_constraint:
            raise ValueError('Dynamic flux value should be set first!')

    def __copy__(self):
        return BaseSampler(
            variable_num=self.variable_num, eq_matrix=self.eq_matrix.copy(),
            eq_right_side_vector=self.eq_right_side_vector.copy(),
            smaller_eq_matrix=self.smaller_eq_matrix.copy(),
            smaller_eq_right_side_vector=self.smaller_eq_right_side_vector.copy(),
            min_value_vector=self.min_value_vector.copy(),
            max_value_vector=self.max_value_vector.copy(),
            dynamic_constraint_index_dict=dict(self.dynamic_constraint_index_dict),
            numeric_eps=self.numeric_eps, verbose=self.verbose)


class LPSampler(BaseSampler):
    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, random_min_factor=default_random_min_factor,
            random_max_factor=default_random_max_factor, maximal_failed_time=default_maximal_failed_time,
            dynamic_constraint_index_dict=None, numeric_eps=eps, verbose=False, *, solver_obj=None):
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        if solver_obj is not None:
            if solver_obj.complete_right_side_array is None:
                raise ValueError()
            super(LPSampler, self).__init__(
                solver_obj.variable_num, solver_obj.complete_flux_constraint_matrix,
                solver_obj.complete_right_side_array, smaller_eq_matrix=None, smaller_eq_right_side_vector=None,
                min_value_vector=solver_obj.min_bound_vector, max_value_vector=solver_obj.max_bound_vector,
                dynamic_constraint_index_dict=dynamic_constraint_index_dict,
                numeric_eps=None, verbose=solver_obj.verbose)
        elif variable_num is not None:
            super(LPSampler, self).__init__(
                variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
                min_value_vector, max_value_vector, dynamic_constraint_index_dict, numeric_eps, verbose)
        else:
            raise ValueError()
        self.random_min_factor = random_min_factor
        self.random_max_factor = random_max_factor
        self.maximal_failed_time = maximal_failed_time

    def refresh(self):
        """
        No need to do anything
        """
        pass

    def set_random_factor(self, new_random_min_factor=None, new_random_max_factor=None):
        if new_random_min_factor is not None:
            self.random_min_factor = new_random_min_factor
        if new_random_max_factor is not None:
            self.random_max_factor = new_random_max_factor

    def set_right_side_vector(self, dynamic_flux_value_dict):
        super(LPSampler, self).set_right_side_vector(dynamic_flux_value_dict)
        # if eq_or_smaller_eq == ParameterName.eq_right_side:
        #     target_vector = self.eq_right_side_vector
        # elif eq_or_smaller_eq == ParameterName.smaller_or_eq_right_side:
        #     target_vector = self.smaller_eq_right_side_vector
        # else:
        #     raise ValueError()
        target_vector = self.eq_right_side_vector
        for flux_name, flux_value in dynamic_flux_value_dict.items():
            target_vector[self.dynamic_constraint_index_dict[flux_name]] = flux_value

    def _clip_to_bound(self, raw_vector):
        clipped_vector = np.clip(raw_vector, self.min_value_vector, self.max_value_vector)
        return clipped_vector

    def _generate_min_factor_bound(self):
        return 1 + np.random.random(self.variable_num) * self.random_min_factor

    def _generate_max_factor_bound(self):
        return np.random.random(self.variable_num) * self.random_max_factor + 1 - self.random_max_factor

    def _generate_lb_ub_for_lp(self):
        # This may generate wrong solution when self.min_value_vector is smaller than 0
        raw_lb_generated_by_min_factor = self.min_value_vector * self._generate_min_factor_bound()
        raw_lb_generated_by_max_factor = self.min_value_vector * self._generate_max_factor_bound()
        raw_lb = np.where(self.min_value_vector >= 0, raw_lb_generated_by_min_factor, raw_lb_generated_by_max_factor)
        clipped_lb = self._clip_to_bound(raw_lb)
        raw_ub = self.max_value_vector * self._generate_max_factor_bound()
        clipped_ub = self._clip_to_bound(raw_ub)
        return clipped_lb, clipped_ub

    def sample(self, n=1):
        super(LPSampler, self).sample(n)
        result_list = []
        failed_time = 0
        while failed_time < n * self.maximal_failed_time:
            lp_lb, lp_ub = self._generate_lb_ub_for_lp()
            if np.all(lp_ub >= lp_lb):
                bounds = np.array([lp_lb, lp_ub]).T
                random_obj = np.random.random(self.eq_matrix.shape[1]) - 0.4
                res = optimize.linprog(
                    random_obj, A_eq=self.eq_matrix, b_eq=self.eq_right_side_vector,
                    A_ub=self.smaller_eq_matrix, b_ub=self.smaller_eq_right_side_vector,
                    bounds=bounds, method="interior-point",
                    options={'tol': self.numeric_eps, })  # "disp": self.verbose
                if res.success:
                    result_list.append(np.array(res.x))
                    if len(result_list) == n:
                        break
            failed_time += 1
        if len(result_list) < n:
            return None
        elif n == 1:
            return result_list[0]
        else:
            return result_list

    def parallel_sample(self, n, processes_num=6, parallel_test=False, display_progress_bar=False, name=None):
        return lp_parallel_sample(
            n, processes_num=processes_num, variable_num=self.variable_num,
            eq_matrix=self.eq_matrix, eq_right_side_vector=self.eq_right_side_vector,
            smaller_eq_matrix=self.smaller_eq_matrix, smaller_eq_right_side_vector=self.smaller_eq_right_side_vector,
            min_value_vector=self.min_value_vector, max_value_vector=self.max_value_vector,
            random_min_factor=self.random_min_factor, random_max_factor=self.random_max_factor,
            numeric_eps=self.numeric_eps, parallel_test=parallel_test, display_progress_bar=display_progress_bar,
            name=name)


def lp_single_sample(parameter_list):
    (
        total_bounds_obj_num, target_num, all_bounds_array, all_obj_array,
        eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
        numeric_eps, send_pipe) = parameter_list
    result_list = []
    for i in range(total_bounds_obj_num):
        res = optimize.linprog(
            all_obj_array[i], A_eq=eq_matrix, b_eq=eq_right_side_vector,
            A_ub=smaller_eq_matrix, b_ub=smaller_eq_right_side_vector,
            bounds=all_bounds_array[i], method="interior-point",
            options={'tol': numeric_eps, })  # "disp": self.verbose
        if res.success:
            result_list.append(np.array(res.x))
            if send_pipe is not None:
                send_pipe.send(1)
            if len(result_list) == target_num:
                break
    return result_list


def generate_valid_lb_ub_bounds_and_obj(
        total_bound_obj_num, variable_num, min_value_vector, max_value_vector, random_min_factor, random_max_factor,
        random_seed):

    def generate_min_factor_bound():
        return 1 + random_seed.random([total_bound_obj_num, variable_num]) * random_min_factor

    def generate_max_factor_bound():
        return random_seed.random([total_bound_obj_num, variable_num]) * random_max_factor + (1 - random_max_factor)

    def clip_to_bound(raw_vector):
        clipped_vector = np.clip(raw_vector, min_value_vector, max_value_vector)
        return clipped_vector

    # This may generate wrong solution when self.min_value_vector is smaller than 0
    raw_lb_generated_by_min_factor = min_value_vector * generate_min_factor_bound()
    raw_lb_generated_by_max_factor = min_value_vector * generate_max_factor_bound()
    raw_lb = np.where(min_value_vector >= 0, raw_lb_generated_by_min_factor, raw_lb_generated_by_max_factor)
    clipped_lb = clip_to_bound(raw_lb)
    raw_ub = max_value_vector * generate_max_factor_bound()
    clipped_ub = clip_to_bound(raw_ub)
    bounds = np.transpose(np.array([clipped_lb, clipped_ub]), axes=(1, 2, 0))
    obj_array = random_seed.random([total_bound_obj_num, variable_num]) - 0.4
    return bounds, obj_array


def lp_parallel_sample(
        target_size, processes_num, variable_num,
        eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,

        min_value_vector, max_value_vector, random_min_factor, random_max_factor, numeric_eps, parallel_test,
        display_progress_bar=False, name=None):

    def parameter_list_generator(current_processes_num, send_pipe):
        for i in range(current_processes_num):
            this_process_target_size = each_process_size_list[i]
            bound_obj_start = i * each_process_valid_bound_num
            bound_obj_end = (i + 1) * each_process_valid_bound_num
            current_bound_array = bounds[bound_obj_start:bound_obj_end]
            current_obj_array = obj_array[bound_obj_start:bound_obj_end]
            yield (
                each_process_valid_bound_num, this_process_target_size, current_bound_array, current_obj_array,
                eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
                numeric_eps, send_pipe)

    def process_result(current_raw_result):
        folded_initial_vector_list.append(current_raw_result)
        return len(current_raw_result)

    # random_seed = np.random.default_rng(4536251)
    random_seed = np.random.default_rng()

    folded_initial_vector_list = []
    current_target_size = target_size
    each_process_minimal_num = 10
    max_iter_num = 100
    iter_count = 0

    if name is not None:
        display_title = f'Initial solution generation for {name}'
    else:
        display_title = 'Initial solution generation'

    with ProgressBarExtraProcess(
        display_progress_bar=display_progress_bar, target_size=target_size, display_title=display_title
    ) as progress_bar:

        # print(f'Computation started')
        while current_target_size > 0 and iter_count < max_iter_num:
            if current_target_size <= each_process_minimal_num:
                parallel_test = True
                current_processes_num = 1
            elif current_target_size <= each_process_minimal_num * processes_num:
                current_processes_num = current_target_size // each_process_minimal_num
            else:
                current_processes_num = processes_num

            each_process_size_list = split_total_num_to_process(current_target_size, current_processes_num)
            finished_num = 0
            each_process_valid_bound_num = each_process_size_list[0] * 2
            valid_bound_num = each_process_valid_bound_num * current_processes_num
            bounds, obj_array = generate_valid_lb_ub_bounds_and_obj(
                valid_bound_num, variable_num, min_value_vector, max_value_vector, random_min_factor, random_max_factor,
                random_seed)
            parameter_list_iter = parameter_list_generator(current_processes_num, progress_bar.send_pipe)

            # print(f'Computation in cycle {iter_count} started')
            if parallel_test:
                for parameter_list in parameter_list_iter:
                    raw_result = lp_single_sample(parameter_list)
                    finished_num += process_result(raw_result)
            else:
                with mp.Pool(processes=current_processes_num) as pool:
                    raw_result_iter = pool.imap(lp_single_sample, parameter_list_iter)
                    for raw_result in raw_result_iter:
                        finished_num += process_result(raw_result)

            current_target_size -= finished_num
            iter_count += 1
            # print(f'Cycle {iter_count} finished')

    # print(f'Sample finished')
    # exit(0)

    if current_target_size > 0:
        return None
    else:
        final_initial_array = np.concatenate(folded_initial_vector_list)[:target_size]
        return final_initial_array
