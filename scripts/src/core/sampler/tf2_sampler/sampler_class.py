from ...common.config import ParamName
from ...common.packages import np, warnings
from ...common.tf_packages import tf, tensor_obj, int_type, float_type, device, tf_function, tf_function_mode

from ..base_class import default_random_min_factor, default_random_max_factor, \
    default_maximal_failed_time, eps, LPSampler
from .common_sampler_functions import pick_random_initial_points, random_direction_pool_generator, \
    random_center_direction_generator, update_saved_tensor, new_step_function, sample_function


class TF2Sampler(object):
    def __init__(
            self, variable_num, eq_matrix=None, eq_right_side_vector=None,
            smaller_eq_matrix=None, smaller_eq_right_side_vector=None,
            min_value_vector=None, max_value_vector=None, dynamic_constraint_index_dict=None,
            numeric_eps=eps, verbose=False):
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        with device:
            if eq_matrix is not None:
                eq_matrix = tensor_obj.constant(eq_matrix)
            if eq_right_side_vector is not None:
                eq_right_side_vector = tensor_obj.variable(eq_right_side_vector, trainable=False)
            if smaller_eq_matrix is not None:
                smaller_eq_matrix = tensor_obj.constant(smaller_eq_matrix)
            if smaller_eq_right_side_vector is not None:
                smaller_eq_right_side_vector = tensor_obj.constant(smaller_eq_right_side_vector)
            if min_value_vector is not None:
                min_value_vector = tensor_obj.constant(min_value_vector)
            if max_value_vector is not None:
                max_value_vector = tensor_obj.constant(max_value_vector)

            self.variable_num = variable_num
            self.eq_matrix = eq_matrix
            self.eq_right_side_vector = eq_right_side_vector
            self.smaller_eq_matrix = smaller_eq_matrix
            self.smaller_eq_right_side_vector = smaller_eq_right_side_vector
            self.min_value_vector = min_value_vector
            self.max_value_vector = max_value_vector
            self.numeric_eps = numeric_eps
            self.verbose = verbose
            self.dynamic_constraint_index_dict = dynamic_constraint_index_dict
            self._set_dynamic_constraint = len(dynamic_constraint_index_dict) == 0

    def set_right_side_vector(self, dynamic_flux_value_dict):
        self._set_dynamic_constraint = True
        with device:
            target_vector = self.eq_right_side_vector
            current_index_list = []
            flux_value_list = []
            for flux_name, flux_value in dynamic_flux_value_dict.items():
                current_index_list.append(self.dynamic_constraint_index_dict[flux_name])
                flux_value_list.append(flux_value)
            target_vector.scatter_nd_update(current_index_list, flux_value_list)

    def sample(self, n=None):
        if not self._set_dynamic_constraint:
            raise ValueError('Dynamic flux value should be set first!')


class TF2HRSampler(TF2Sampler):
    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, sampler_options_dict,
            dynamic_constraint_index_dict=None, verbose=False, refresh=True):
        if dynamic_constraint_index_dict is None:
            dynamic_constraint_index_dict = {}
        with device:
            initial_random_min_factor = sampler_options_dict.get_option(
                ParamName.initial_random_min_factor, default_random_min_factor)
            initial_random_max_factor = sampler_options_dict.get_option(
                ParamName.initial_random_max_factor, default_random_max_factor)
            initial_maximal_failed_time = sampler_options_dict.get_option(
                ParamName.initial_maximal_failed_time, default_maximal_failed_time)
            numeric_eps = sampler_options_dict.get_option(ParamName.numeric_eps, eps)

            super(TF2HRSampler, self).__init__(
                variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
                min_value_vector, max_value_vector, dynamic_constraint_index_dict, numeric_eps, verbose)

            self.lp_sampler = LPSampler(
                variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
                min_value_vector, max_value_vector, random_min_factor=initial_random_min_factor,
                random_max_factor=initial_random_max_factor, maximal_failed_time=initial_maximal_failed_time,
                dynamic_constraint_index_dict=dynamic_constraint_index_dict, numeric_eps=numeric_eps, verbose=verbose)
            self.thinning = sampler_options_dict.get_option(ParamName.thinning, 50)
            self.batch_size = sampler_options_dict.get_option(ParamName.batch_size, 50)

            self.projection_tensor = tensor_obj.constant(projection_matrix)

            self._warmup_tensor = None
            warmup_multiplier = sampler_options_dict.get_option(ParamName.warmup_multiplier, 4)
            self._total_lp_warmup_num = warmup_multiplier * self.batch_size
            self._warmup_mix_num = warmup_multiplier

            self._random_direction_pool = None
            random_direction_pool_multiplier = sampler_options_dict.get_option(
                ParamName.random_direction_pool_multiplier, 100)
            self._random_direction_pool_size = random_direction_pool_multiplier * self.batch_size
            self._maximal_direction_sample_num = self._random_direction_pool_size / 2
            self._current_direction_sample_num = None

            # For numerical stability
            self.maximal_range_ratio = sampler_options_dict.get_option(ParamName.maximal_range_ratio, 0.95)
            self._maximal_fail_num = sampler_options_dict.get_option(ParamName.maximal_fail_num, 10)
            # self._lp_warmup()
            if refresh:
                TF2HRSampler.refresh(self)
                # self._update_random_pool()

    def _update_random_pool(self):
        self._random_direction_pool = random_direction_pool_generator(
            self.variable_num, self._random_direction_pool_size, self.projection_tensor, self.numeric_eps)
        self._current_direction_sample_num = 0

    def refresh(self):
        with device:
            self._lp_warmup()
            self._update_random_pool()
            pass

    def _lp_warmup(self):
        with device:
            failed_count = 0
            warmup_list = self.lp_sampler.sample(self._total_lp_warmup_num)
            while warmup_list is None and failed_count < 3:
                failed_count += 1
                warnings.warn("Cannot generate enough warm-up solutions. Relax random condition")
                self.lp_sampler.set_random_factor(
                    new_random_min_factor=self.lp_sampler.random_min_factor / 2,
                    new_random_max_factor=self.lp_sampler.random_max_factor / 2)
                warmup_list = self.lp_sampler.sample(self._total_lp_warmup_num)
            if warmup_list is None:
                raise ValueError("Cannot generate enough warm-up solutions")
            self._warmup_tensor = tensor_obj.constant(warmup_list)

    def random_warmup_point(self, points_num=None):
        with device:
            if points_num is None:
                points_num = self.batch_size
            return pick_random_initial_points(
                points_num, self._total_lp_warmup_num, self._warmup_tensor, self._warmup_mix_num)

    def generate_random_direction(self, direction_num=None):
        with device:
            if direction_num is None:
                direction_num = self.batch_size
            if self._current_direction_sample_num >= self._maximal_direction_sample_num:
                self._update_random_pool()
            self._current_direction_sample_num += direction_num
            random_col_index_array = tensor_obj.random_uniform(
                (direction_num,), dtype=int_type, maxval=self._random_direction_pool_size)
            target_random_direction_tensor = tf.gather(self._random_direction_pool, random_col_index_array)
            return target_random_direction_tensor

    def set_right_side_vector(self, dynamic_flux_value_dict):
        super(TF2HRSampler, self).set_right_side_vector(dynamic_flux_value_dict)
        self.lp_sampler.set_right_side_vector(dynamic_flux_value_dict)
        self.refresh()

    def sample(self, n=None):
        super(TF2HRSampler, self).sample(n)


class TF2OptGpSampler(TF2HRSampler):
    def __init__(
            self, variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
            min_value_vector, max_value_vector, projection_matrix, sampler_options_dict,
            dynamic_constraint_index_dict=None, verbose=False):
        with device:
            super(TF2OptGpSampler, self).__init__(
                variable_num, eq_matrix, eq_right_side_vector, smaller_eq_matrix, smaller_eq_right_side_vector,
                min_value_vector, max_value_vector, projection_matrix, sampler_options_dict,
                dynamic_constraint_index_dict=dynamic_constraint_index_dict, verbose=verbose, refresh=False)

            self._maximal_warmup_steps = sampler_options_dict.get_option('maximal_warmup_steps', 100)
            self._maximal_sample_index = sampler_options_dict.get_option('maximal_sample_index', 10000)
            saved_sample_num_multiplier = sampler_options_dict.get_option('saved_sample_num_multiplier', 50)
            self._saved_sample_num = saved_sample_num_multiplier * self.batch_size
            self._saved_sample_update_rate = sampler_options_dict.get_option('saved_sample_update_rate', 0.05)

            self._sum_num_for_center = 0
            self.zero_batch_size_np_array = np.zeros(self.batch_size)
            self._current_sum_tensor = tensor_obj.variable(
                np.zeros(self.variable_num), trainable=False)
            self._current_center = tensor_obj.variable(
                np.zeros(self.variable_num), trainable=False)
            self._current_location_tensor = tensor_obj.variable(
                np.zeros((self.batch_size, self.variable_num)), trainable=False)
            self._saved_sample_tensor = tensor_obj.variable(
                np.zeros((self._saved_sample_num, self.variable_num)), trainable=False)
            self._failed_count_tensor = tensor_obj.variable(
                self.zero_batch_size_np_array, dtype=int_type, trainable=False)
            self._step_count_tensor = tensor_obj.variable(
                self.zero_batch_size_np_array, dtype=int_type, trainable=False)
            TF2OptGpSampler.refresh(self)

    def refresh(self):
        with device:
            super().refresh()
            new_batch_tensor = self.random_warmup_point()
            self._current_location_tensor.assign(new_batch_tensor)
            self._sum_num_for_center = self.batch_size
            self._current_sum_tensor.assign(tf.reduce_sum(new_batch_tensor, axis=0))
            self._current_center.assign(self._current_sum_tensor / self._sum_num_for_center)
            self._saved_sample_tensor.assign(
                tf.gather(new_batch_tensor, np.random.choice(self.batch_size, size=self._saved_sample_num)))
            self._failed_count_tensor.assign(self.zero_batch_size_np_array)
            self._step_count_tensor.assign(self.zero_batch_size_np_array)

    @tf_function(tf_function_mode)
    def _generate_random_center_direction(self, direction_num=None):
        with device:
            if direction_num is None:
                direction_num = self.batch_size
            return random_center_direction_generator(
                self._saved_sample_tensor, self._saved_sample_num, self._current_center,
                direction_num, self.projection_tensor, eps)

    @tf_function(tf_function_mode)
    def _update_center(self, new_sampled_points, new_point_size):
        with device:
            self._sum_num_for_center += new_point_size
            self._current_sum_tensor.assign_add(tf.reduce_sum(new_sampled_points, axis=0))
            self._current_center.assign(self._current_sum_tensor / tf.cast(self._sum_num_for_center, dtype=float_type))

    @tf_function(tf_function_mode)
    def _update_saved_tensor(self, new_sampled_points, new_point_size):
        with device:
            update_saved_tensor(
                self._saved_sample_tensor, self._saved_sample_num, self._saved_sample_update_rate,
                new_point_size, new_sampled_points)

    def new_step(self, pickout_tensor_storage_func):
        with device:
            new_step_function(
                self, self._current_location_tensor, self._step_count_tensor, self._failed_count_tensor,
                self._maximal_warmup_steps, self._maximal_sample_index, self._maximal_fail_num, self.thinning,
                self.batch_size, pickout_tensor_storage_func)

    def sample(self, n=None):
        super(TF2HRSampler, self).sample(n)
        with device:
            if n is None:
                n = 2 * self.batch_size
            output_sample_tensor = tensor_obj.variable(np.zeros((n, self.variable_num)), trainable=False)
            current_size_tensor = tensor_obj.variable(0, dtype=int_type, trainable=False)
            sample_function(self, n, output_sample_tensor, current_size_tensor)
            return output_sample_tensor
