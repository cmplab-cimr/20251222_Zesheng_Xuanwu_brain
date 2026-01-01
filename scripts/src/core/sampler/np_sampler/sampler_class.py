from ...common.classes import CircularQueue
from ...common.packages import np, warnings, mp
from ...common.config import ParamName
from ...common.functions import split_total_num_to_process, ProgressBarExtraProcess

from ..base_class import default_random_min_factor, default_random_max_factor, \
    default_maximal_failed_time, eps, BaseSampler, LPSampler
from .common_sampler_functions import calculate_direction_range, check_new_location, \
    random_direction_pool_generator

rand_rng = np.random.default_rng(4536251)


class ParallelHRSampler(BaseSampler):
    """The abstract base class for hit-and-run samplers.

    Parameters
    ----------
    thinning : int
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.

    Attributes
    ----------
    thinning : int
        The currently used thinning factor.

    """

    def __init__(self, base_lp_sampler, thinning=50, processes_num=10, batch_size=20):
        super().__init__(
            base_lp_sampler.variable_num, base_lp_sampler.eq_matrix, base_lp_sampler.eq_right_side_vector,
            base_lp_sampler.smaller_eq_matrix, base_lp_sampler.smaller_eq_right_side_vector,
            base_lp_sampler.min_value_vector, base_lp_sampler.max_value_vector,
            dynamic_constraint_index_dict=base_lp_sampler.dynamic_constraint_index_dict,
            numeric_eps=base_lp_sampler.numeric_eps, verbose=base_lp_sampler.verbose)

        self.lp_sampler = base_lp_sampler

        self.valid_smaller_eq_matrix = None
        self.valid_transformed_smaller_eq_matrix = None
        self.valid_smaller_eq_right_side_vector = None
        self.valid_flux_constraint_matrix_q = None
        self.free_flux_num = None

        self.thinning = thinning
        self.processes_num = processes_num
        self.batch_size = batch_size
        self.normal_a_min = 0.5
        self.normal_a_max = 1.5

    def initialize(self):
        (
            self.valid_smaller_eq_matrix, self.valid_transformed_smaller_eq_matrix,
            self.valid_smaller_eq_right_side_vector,
            self.valid_flux_constraint_matrix_q, self.free_flux_num
        ) = initialize_equality_inequality_constraints(
            self.eq_matrix, self.eq_right_side_vector,
            self.smaller_eq_matrix, self.smaller_eq_right_side_vector,
            self.min_value_vector, self.max_value_vector, self.numeric_eps)

    def sample(self, n=1, parallel_test=False, display_progress_bar=True):
        return parallel_sample(
            self.lp_sampler, n,
            self.valid_smaller_eq_matrix, self.valid_transformed_smaller_eq_matrix,
            self.valid_smaller_eq_right_side_vector,
            self.valid_flux_constraint_matrix_q, self.processes_num,
            self.batch_size, self.free_flux_num, self.thinning,
            self.numeric_eps, self.normal_a_min, self.normal_a_max,
            parallel_test=parallel_test, display_progress_bar=display_progress_bar,
        )


def initialize_equality_inequality_constraints(
        complete_flux_constraint_matrix, complete_right_side_array,
        smaller_eq_matrix, smaller_eq_right_side_vector, min_value_vector, max_value_vector, computation_eps):
    """
    A @ (x + x0) = b
    C @ (x + x0) <= d

    A @ x = 0
    C @ x <= d - C @ x0 = d'

    LQ decomposition to A: A = L @ Q

    L @ Q @ x = 0
    y = Q @ x, x = Q.T @ y
    L @ y = 0

    C @ Q.T @ y <= d'
    C' @ y <= d'

    """
    from ...solver.solver_construction_functions.common_construct_functions import (
        constraint_matrix_verification_and_simplification)

    (
        valid_constraint_matrix, valid_right_side_array, valid_constraint_num
    ) = constraint_matrix_verification_and_simplification(
        complete_flux_constraint_matrix, complete_right_side_array, computation_eps)

    variable_num = complete_flux_constraint_matrix.shape[1]
    flux_constraint_matrix_q_t, matrix_r = np.linalg.qr(valid_constraint_matrix.T, mode='complete')
    valid_flux_constraint_matrix_q_t = flux_constraint_matrix_q_t[:, valid_constraint_num:]
    free_flux_num = variable_num - valid_constraint_num
    # transformed_smaller_eq_matrix = smaller_eq_matrix @ valid_flux_constraint_matrix_q_t
    variable_identity_matrix = np.identity(variable_num, dtype='float64')
    complete_smaller_eq_matrix_list = [
        variable_identity_matrix,
        -variable_identity_matrix
    ]
    if smaller_eq_matrix is not None:
        complete_smaller_eq_matrix_list.append(smaller_eq_matrix)
    complete_smaller_eq_right_side_list = [
        max_value_vector,
        -min_value_vector
    ]
    if smaller_eq_right_side_vector is not None:
        complete_smaller_eq_right_side_list.append(smaller_eq_right_side_vector)
    valid_flux_constraint_matrix_q = valid_flux_constraint_matrix_q_t.T
    complete_smaller_eq_matrix = np.vstack(complete_smaller_eq_matrix_list)
    complete_smaller_eq_right_side_vector = np.concatenate(complete_smaller_eq_right_side_list)
    transformed_smaller_eq_matrix = (complete_smaller_eq_matrix @ valid_flux_constraint_matrix_q_t)
    # (
    #     valid_smaller_eq_matrix, valid_transformed_smaller_eq_matrix, valid_smaller_eq_right_side_array
    # ) = smaller_eq_bound_matrix_simplification(
    #     complete_smaller_eq_matrix, valid_flux_constraint_matrix_q_t, complete_smaller_eq_right_side_vector, computation_eps)
    return (
        complete_smaller_eq_matrix, transformed_smaller_eq_matrix, complete_smaller_eq_right_side_vector,
        valid_flux_constraint_matrix_q, free_flux_num)


def single_vector_sample(parameter_list):

    def find_max_ratio(y_matrix, p_matrix, matrix_c_t, right_side_matrix_d):
        """
        C' @ (y + a * p) <= d'

        a * (C @ p) <= d - C @ y
        """
        delta_d = right_side_matrix_d - y_matrix @ matrix_c_t
        projected_p = p_matrix @ matrix_c_t
        raw_a_matrix = delta_d / (projected_p + 0.532718273 * computation_eps)
        raw_a_matrix[raw_a_matrix < 0] = 9e10
        raw_a_matrix[raw_a_matrix < computation_eps] = computation_eps
        each_max_a_max = np.min(raw_a_matrix, axis=1, keepdims=True)
        return each_max_a_max

    (
        complete_initial_solution_array, complete_smaller_eq_right_side_array,
        target_size, transformed_smaller_eq_matrix_t, smaller_eq_right_side_vector,
        valid_flux_constraint_matrix_q,
        batch_size, free_flux_num, thinning, normal_a_min, normal_a_max, computation_eps, send_pipe
    ) = parameter_list

    max_continuous_invalid_step = 20
    max_valid_step = 100000
    minimal_step_norm = 1e-5
    max_a_ratio = 0.8
    valid_step_count = np.zeros(batch_size, dtype=int)
    continuous_invalid_step_count = np.zeros(batch_size, dtype=int)
    selected_solution = np.zeros(batch_size, dtype=bool)
    output_free_flux_list = []
    output_initial_flux_list = []
    complete_initial_solution_size = complete_initial_solution_array.shape[0]
    rest_initial_num = complete_initial_solution_size

    current_free_flux_array = np.zeros([batch_size, free_flux_num], dtype='float64')
    current_initial_flux_array = complete_initial_solution_array[:batch_size]
    current_smaller_eq_right_side_array = complete_smaller_eq_right_side_array[:batch_size]
    random_free_direction_zero_matrix = np.zeros(free_flux_num, dtype ='float64')
    random_free_direction_eye_matrix = np.eye(free_flux_num, dtype ='float64')

    while len(output_free_flux_list) < target_size:
        random_free_direction_vector = rand_rng.multivariate_normal(
            random_free_direction_zero_matrix, random_free_direction_eye_matrix, batch_size)
        random_free_direction_vector /= np.linalg.norm(random_free_direction_vector, axis=1, keepdims=True)

        each_free_flux_max_a_max = find_max_ratio(
            current_free_flux_array, random_free_direction_vector,
            transformed_smaller_eq_matrix_t, current_smaller_eq_right_side_array)
        real_each_free_flux_max_a_max = each_free_flux_max_a_max * max_a_ratio

        each_free_flux_a_min = np.minimum(normal_a_min, real_each_free_flux_max_a_max)
        each_free_flux_a_max = np.minimum(normal_a_max, real_each_free_flux_max_a_max)
        each_free_flux_a = rand_rng.random([batch_size, 1]) * (each_free_flux_a_max - each_free_flux_a_min) + each_free_flux_a_min

        new_free_flux_array = current_free_flux_array + each_free_flux_a * random_free_direction_vector
        delta_norm = np.linalg.norm((new_free_flux_array - current_free_flux_array), axis=1)
        current_step_valid_array = delta_norm > minimal_step_norm

        valid_step_count[current_step_valid_array] += 1
        continuous_invalid_step_count[current_step_valid_array] = 0
        continuous_invalid_step_count[~current_step_valid_array] += 1
        selected_solution[current_step_valid_array] = False

        success_array = (
            (valid_step_count > 0)
            & ((valid_step_count % thinning) == 0)
            & ~selected_solution)
        if np.any(success_array):
            success_free_flux = current_free_flux_array[success_array]
            success_size = success_free_flux.shape[0]
            output_free_flux_list.extend(success_free_flux)
            output_initial_flux_list.extend(current_initial_flux_array[success_array])
            selected_solution[success_array] = True
            if send_pipe is not None:
                send_pipe.send(success_size)
        failed_array = (
            (continuous_invalid_step_count > max_continuous_invalid_step)
            | (valid_step_count > max_valid_step))
        if np.any(failed_array):
            failed_num = np.count_nonzero(failed_array)
            if rest_initial_num < failed_num:
                rest_failed_num = failed_num - rest_initial_num
                new_initial_free_flux_array = np.vstack([
                    complete_initial_solution_array[-rest_initial_num:],
                    complete_initial_solution_array[:rest_failed_num]
                ])
                new_smaller_eq_right_side_array = np.vstack([
                    complete_smaller_eq_right_side_array[-rest_initial_num:],
                    complete_smaller_eq_right_side_array[:rest_failed_num]
                ])
                rest_initial_num = complete_initial_solution_size - rest_failed_num
            else:
                initial_start = -rest_initial_num
                initial_end = initial_start + failed_num
                new_initial_free_flux_array = complete_initial_solution_array[initial_start:initial_end]
                new_smaller_eq_right_side_array = complete_smaller_eq_right_side_array[initial_start:initial_end]
                rest_initial_num -= failed_num
            current_free_flux_array[failed_array] = new_initial_free_flux_array
            current_initial_flux_array[failed_array] = new_initial_free_flux_array
            current_smaller_eq_right_side_array[failed_array] = new_smaller_eq_right_side_array
            valid_step_count[failed_array] = 0
            continuous_invalid_step_count[failed_array] = 0

    output_initial_flux_array = np.array(output_initial_flux_list)
    output_free_flux_array = np.array(output_free_flux_list)
    output_original_flux_array = output_initial_flux_array + output_free_flux_array @ valid_flux_constraint_matrix_q
    return output_original_flux_array


def parallel_sample(
        lp_sampler, target_size, smaller_eq_matrix, transformed_smaller_eq_matrix, smaller_eq_right_side_vector,
        valid_flux_constraint_matrix_q, processes_num,
        batch_size, free_flux_num, thinning, computation_eps, normal_a_min, normal_a_max,
        parallel_test=False, display_progress_bar=True, average_initial_num=3
):
    def generate_new_initial(required_initial_num):
        """Row-based solution"""
        raw_initial_solution_array = lp_sampler.parallel_sample(
            required_initial_num, processes_num=processes_num, parallel_test=parallel_test)
        random_index = rand_rng.choice(required_initial_num, size=(required_initial_num, average_initial_num))
        new_initial_solution_array = raw_initial_solution_array[random_index].mean(axis=1)
        new_smaller_eq_right_side_array = (
            np.reshape(smaller_eq_right_side_vector, [1, -1]) -
            (smaller_eq_matrix @ new_initial_solution_array.T).T)
        return new_initial_solution_array, new_smaller_eq_right_side_array

    def process_result(current_raw_result):
        folded_flux_array_list.append(current_raw_result)

    def parameter_list_generator(_target_size, _batch_size, _process_num):
        each_initial_num = _batch_size * 2
        complete_initial_num = each_initial_num * _process_num
        transformed_smaller_eq_matrix_t = transformed_smaller_eq_matrix.T
        each_process_size_list = split_total_num_to_process(_target_size, _process_num)
        (
            complete_initial_solution_array, complete_smaller_eq_right_side_array
        ) = generate_new_initial(complete_initial_num)
        for i in range(_process_num):
            this_process_target_size = each_process_size_list[i]
            bound_obj_start = i * each_initial_num
            bound_obj_end = (i + 1) * each_initial_num
            current_initial_solution_array = complete_initial_solution_array[bound_obj_start:bound_obj_end]
            current_smaller_eq_right_side_array = complete_smaller_eq_right_side_array[bound_obj_start:bound_obj_end]
            yield (
                current_initial_solution_array, current_smaller_eq_right_side_array,
                this_process_target_size, transformed_smaller_eq_matrix_t, smaller_eq_right_side_vector,
                valid_flux_constraint_matrix_q,
                batch_size, free_flux_num, thinning, normal_a_min, normal_a_max, computation_eps, send_pipe)

    folded_flux_array_list = []
    with ProgressBarExtraProcess(
            display_progress_bar=display_progress_bar, target_size=target_size, display_title='Sampling'
    ) as progress_bar:
        send_pipe = progress_bar.send_pipe
        parameter_list_iter = parameter_list_generator(target_size, batch_size, processes_num)
        if parallel_test:
            for parameter_list in parameter_list_iter:
                raw_result = single_vector_sample(parameter_list)
                process_result(raw_result)
        else:
            with mp.Pool(processes=processes_num) as pool:
                raw_result_iter = pool.imap(single_vector_sample, parameter_list_iter)
                for raw_result in raw_result_iter:
                    process_result(raw_result)

    final_flux_solution_array = np.concatenate(folded_flux_array_list)
    return final_flux_solution_array


class HRSampler(BaseSampler):
    """The abstract base class for hit-and-run samplers.

    Parameters
    ----------
    thinning : int
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.

    Attributes
    ----------
    thinning : int
        The currently used thinning factor.

    """

    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, thinning=50,
            initial_random_min_factor=default_random_min_factor, initial_random_max_factor=default_random_max_factor,
            initial_maximal_failed_time=default_maximal_failed_time,
            dynamic_constraint_index_dict=None, numeric_eps=eps, verbose=False):
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        super(HRSampler, self).__init__(
            variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, dynamic_constraint_index_dict=dynamic_constraint_index_dict,
            numeric_eps=numeric_eps, verbose=verbose)

        self.lp_sampler = LPSampler(
            variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, random_min_factor=initial_random_min_factor,
            random_max_factor=initial_random_max_factor, maximal_failed_time=initial_maximal_failed_time,
            dynamic_constraint_index_dict=dynamic_constraint_index_dict, numeric_eps=numeric_eps, verbose=verbose)
        self.thinning = thinning

        self.warmup_array = None  # Need refresh
        self.warmup_num = 100

        self.projection_matrix = projection_matrix

        self.random_direction_pool = None  # Need refresh
        self.random_direction_pool_size = 1000
        self.maximal_direction_sample_num = round(self.random_direction_pool_size * np.sqrt(self.variable_num))
        # self.current_direction_sample_num = self.maximal_direction_sample_num  # Need refresh
        self.current_direction_sample_num = None  # Need refresh

        self.maximal_range_ratio = 0.95  # For numerical stability

        self.maximal_trial_num = 15
        self.recent_points_num = 20
        self.recent_points_store_interval = 100
        # self.recent_points_circular_queue = CircularQueue(self.recent_points_num)
        self.recent_points_circular_queue = None  # Need refresh
        HRSampler.refresh(self)

    def refresh(self):
        self._lp_warmup()
        self.random_direction_pool = None
        self.current_direction_sample_num = self.maximal_direction_sample_num
        self.recent_points_circular_queue = CircularQueue(self.recent_points_num)

    def _lp_parallel_generator(self, total_num, process_num):
        sample_num_list = split_total_num_to_process(
            total_num=total_num, processes_num=process_num)

    def _lp_warmup(self):
        failed_count = 0
        warmup_list = self.lp_sampler.sample(self.warmup_num)
        while warmup_list is None and failed_count < 3:
            failed_count += 1
            warnings.warn("Cannot generate enough warm-up solutions. Relax random condition")
            self.lp_sampler.set_random_factor(
                new_random_min_factor=self.lp_sampler.random_min_factor / 2,
                new_random_max_factor=self.lp_sampler.random_max_factor / 2)
            warmup_list = self.lp_sampler.sample(self.warmup_num)
        if warmup_list is None:
            raise ValueError("Cannot generate enough warm-up solutions")
        self.warmup_array = np.array(warmup_list)

    def random_warmup_point(self):
        """Find an approximately random point in the flux cone."""

        index_array = np.random.randint(
            self.warmup_num, size=min(3, np.ceil(np.sqrt(self.warmup_num)))
        )
        return self.warmup_array[index_array, :].mean(axis=0)

    def generate_random_direction(self):
        if self.current_direction_sample_num >= self.maximal_direction_sample_num:
            self.random_direction_pool = random_direction_pool_generator(
                self.variable_num, self.random_direction_pool_size, self.projection_matrix, self.numeric_eps)
            self.current_direction_sample_num = 0
        self.current_direction_sample_num += 1
        random_col_index = np.random.randint(self.random_direction_pool_size)
        return self.random_direction_pool[:, random_col_index]

    def _store_current_location(self, current_location, sample_index=None):
        pass

    def _trace_back(self):
        pass

    def set_right_side_vector(self, dynamic_flux_value_dict):
        super(HRSampler, self).set_right_side_vector(dynamic_flux_value_dict)
        # if eq_or_smaller_eq == ParameterName.eq_right_side:
        #     target_vector = self.eq_right_side_vector
        # elif eq_or_smaller_eq == ParameterName.smaller_or_eq_right_side:
        #     target_vector = self.smaller_eq_right_side_vector
        # else:
        #     raise ValueError()
        target_vector = self.eq_right_side_vector
        for flux_name, flux_value in dynamic_flux_value_dict.items():
            target_vector[self.dynamic_constraint_index_dict[flux_name]] = flux_value


class OptGpSampler(HRSampler):
    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, thinning=50,
            initial_random_min_factor=default_random_min_factor, initial_random_max_factor=default_random_max_factor,
            initial_maximal_failed_time=default_maximal_failed_time, dynamic_constraint_index_dict=None,
            numeric_eps=eps, verbose=False):
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        super().__init__(
            variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, thinning,
            initial_random_min_factor, initial_random_max_factor, initial_maximal_failed_time,
            dynamic_constraint_index_dict=dynamic_constraint_index_dict, numeric_eps=numeric_eps, verbose=verbose)
        self.current_sample_num = None  # Need refresh
        self.current_sum_vector = None  # Need refresh
        self.current_center = None  # Need refresh
        # self.total_output_sample_list = []  # Need refresh

        self.maximal_sample_index = 10000
        self.saved_sample_num = 200
        self.saved_sample_update_rate = 0.05
        # self.saved_sample_list = list(self.warmup_array)
        self.saved_sample_list = []  # Need refresh
        OptGpSampler.refresh(self)

    def refresh(self):
        super().refresh()
        self.current_sample_num = self.warmup_num
        self.current_sum_vector = np.sum(self.warmup_array, axis=0)
        self.current_center = self.current_sum_vector / self.current_sample_num
        # self.total_output_sample_list = []

        self.saved_sample_list = list(self.warmup_array)

    def _random_center_unit_direction(self):
        raw_direction = self.saved_sample_list[np.random.randint(len(self.saved_sample_list))] - self.current_center
        # Numpy will deal with this vector product automatically
        projected_direction = self.projection_matrix @ raw_direction
        unit_direction = projected_direction / np.linalg.norm(projected_direction)
        return unit_direction

    def _update_center(self, new_sample_point):
        self.current_sample_num += 1
        self.current_sum_vector += new_sample_point
        self.current_center = self.current_sum_vector / self.current_sample_num

    def _trace_back_center(self, added_sample_point):
        self.current_sample_num -= 1
        self.current_sum_vector -= added_sample_point
        self.current_center = self.current_sum_vector / self.current_sample_num

    def _update_saved_sample(self, new_sample_point):
        if len(self.saved_sample_list) < self.saved_sample_num:
            self.saved_sample_list.append(new_sample_point)
        else:
            if np.random.random() < self.saved_sample_update_rate:
                self.saved_sample_list[np.random.randint(len(self.saved_sample_list))] = new_sample_point

    def _new_step(self, current_location):
        final_next_location = None
        # Add a trace back mechanism
        for _ in range(self.maximal_trial_num):
            random_unit_direction = self._random_center_unit_direction()
            direction_range = calculate_direction_range(
                current_location, random_unit_direction, self.min_value_vector, self.max_value_vector,
                self.smaller_eq_matrix, self.smaller_eq_right_side_vector, self.numeric_eps)
            if direction_range is None:
                continue
            random_range = np.random.random() * direction_range * self.maximal_range_ratio
            next_location = current_location + random_unit_direction * random_range
            valid_new_location = check_new_location(
                next_location, self.min_value_vector, self.max_value_vector,
                self.smaller_eq_matrix, self.smaller_eq_right_side_vector, self.numeric_eps)
            if valid_new_location:
                final_next_location = next_location
                break
        return final_next_location

    def _store_current_location(self, current_location, sample_index=None):
        self._update_center(current_location)
        self._update_saved_sample(current_location)
        if sample_index is None or sample_index % self.recent_points_store_interval == 0:
            self.recent_points_circular_queue.push(current_location)

    def _trace_back(self):
        if self.recent_points_circular_queue.is_empty():
            return None
        current_point = self.recent_points_circular_queue.tail_out()
        self._trace_back_center(current_point)
        return current_point

    def sample(self, n=1):
        super(OptGpSampler, self).sample(n)
        output_sample_list = []
        thinning_value = self.thinning
        current_point = self.random_warmup_point()
        sample_index = 1
        while len(output_sample_list) < n:
            if sample_index % thinning_value == 0:
                output_sample_list.append(current_point)
            new_sample_point = self._new_step(current_point)
            if new_sample_point is None:
                current_point = self._trace_back()
            else:
                self._store_current_location(current_point, sample_index)
                current_point = new_sample_point
                sample_index += 1
            if current_point is None or sample_index >= self.maximal_sample_index:
                sample_index = 0
                current_point = self.random_warmup_point()
        return output_sample_list


# Deprecated
class CHRRSampler(HRSampler):
    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, thinning=50,
            initial_random_min_factor=default_random_min_factor, initial_random_max_factor=default_random_max_factor,
            initial_maximal_failed_time=default_maximal_failed_time, numeric_eps=eps, verbose=False):
        super().__init__(
            variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, thinning=thinning,
            initial_random_min_factor=initial_random_min_factor, initial_random_max_factor=initial_random_max_factor,
            initial_maximal_failed_time=initial_maximal_failed_time, numeric_eps=numeric_eps, verbose=verbose)
        pass

