from ...common.packages import np, np_int_type
from ...common.tf_packages import tf, tfp, tensor_obj, int_type, tf_function, tf_function_mode


def generate_random_segment_tensor(batch_size, segment_num):
    separator_num = segment_num - 1
    raw_ordered_random_separator = tf.sort(tensor_obj.random_uniform((batch_size, separator_num)), axis=1)
    right_side_separator = tf.concat([raw_ordered_random_separator, tensor_obj.ones([batch_size, 1])], axis=1)
    left_side_separator = tf.concat([tensor_obj.zeros([batch_size, 1]), raw_ordered_random_separator], axis=1)
    random_segment_tensor = right_side_separator - left_side_separator
    return random_segment_tensor


@tf_function(tf_function_mode)
def pick_random_initial_points(
        batch_size, total_warmup_num, warmup_tensor, warmup_mix_num):
    warmup_index_array = tensor_obj.random_uniform(
        (batch_size, warmup_mix_num, 1), dtype=int_type, maxval=total_warmup_num)
    raw_initial_tensor = tf.gather_nd(warmup_tensor, warmup_index_array)
    warmup_mix_tensor = generate_random_segment_tensor(batch_size, warmup_mix_num)
    mixed_warmup_tensor = tf.reduce_sum(
        raw_initial_tensor * tf.reshape(warmup_mix_tensor, [batch_size, warmup_mix_num, 1]), axis=1)
    return mixed_warmup_tensor


@tf_function(tf_function_mode)
def random_direction_pool_generator(variable_num, pool_size, projection_tensor, eps):
    mean = tensor_obj.zeros(variable_num)
    # cov = tensor_obj.eye(variable_num)
    mvn_obj = tfp.distributions.MultivariateNormalDiag(loc=mean)
    raw_direction_vector_pool = mvn_obj.sample(pool_size)
    if projection_tensor is not None:
        projected_direction_pool = tf.transpose(projection_tensor @ tf.transpose(raw_direction_vector_pool))
    else:
        projected_direction_pool = raw_direction_vector_pool
    normalized_direction_pool = tf.math.l2_normalize(projected_direction_pool, axis=1, epsilon=eps)
    return normalized_direction_pool


@tf_function(tf_function_mode)
def random_center_direction_generator(
        saved_sample_tensor, saved_sample_num, current_center, batch_size, projection_tensor, eps):
    random_indices = np.random.choice(saved_sample_num, batch_size, replace=False)
    raw_direction = tf.gather(saved_sample_tensor, random_indices) - current_center
    if projection_tensor is not None:
        projected_direction = tf.transpose(projection_tensor @ tf.transpose(raw_direction))
    else:
        projected_direction = raw_direction
    normalized_direction = tf.math.l2_normalize(projected_direction, axis=1, epsilon=eps)
    return normalized_direction


@tf_function(tf_function_mode)
def update_saved_tensor(
        saved_sample_tensor, saved_sample_num, saved_sample_update_rate, batch_size, new_sample_batch_point):
    def non_replace_np_choice(_saved_sample_num, _update_size):
        return np.random.choice(_saved_sample_num, size=_update_size, replace=False).astype(np_int_type)

    update_tensor = tensor_obj.random_uniform((batch_size,)) < saved_sample_update_rate
    update_size = tf.math.count_nonzero(update_tensor, dtype=int_type)
    if update_size > tensor_obj.constant(0, dtype=int_type):
        # int_update_size = int(update_size)
        # updated_indices = np.random.choice(saved_sample_num, size=int_update_size, replace=False)
        updated_indices = tf.numpy_function(non_replace_np_choice, [saved_sample_num, update_size], int_type)
        updated_content = tf.boolean_mask(new_sample_batch_point, update_tensor)
        saved_sample_tensor.scatter_nd_update(updated_indices, updated_content)


@tf_function(tf_function_mode)
def generate_random_direction(self, step_count_tensor, maximal_warmup_steps):
    pure_random_unit_direction = self.generate_random_direction()
    center_random_unit_direction = self._generate_random_center_direction()
    final_random_direction = tf.where(
        tf.reshape(step_count_tensor < maximal_warmup_steps, [-1, 1]),
        pure_random_unit_direction, center_random_unit_direction)
    return final_random_direction


@tf_function(tf_function_mode)
def calculate_direction_range(
        batch_size, current_location_tensor, unit_direction, lb, ub, a_ub, b_ub, eps):
    """
    Calculate the feasible maximal range of the given direction from current location.
    The current location should be feasible, the direction should be unit vector and satisfy linear constraints,
        and only positive value will be calculated.
    """
    def filter_and_min(raw_distance_tensor):
        filtered_distance_tensor = tf.where(raw_distance_tensor < 0, tensor_obj.constant(np.inf), raw_distance_tensor)
        min_distance_tensor = tf.reduce_min(filtered_distance_tensor, axis=1)
        return min_distance_tensor

    final_range_tensor = tensor_obj.constant([np.inf] * batch_size)
    flux_lb_distance_tensor = (lb + eps - current_location_tensor) / (unit_direction + eps)
    final_range_tensor = tf.minimum(final_range_tensor, filter_and_min(flux_lb_distance_tensor))

    flux_ub_distance_vector = (ub - eps - current_location_tensor) / (unit_direction + eps)
    final_range_tensor = tf.minimum(final_range_tensor, filter_and_min(flux_ub_distance_vector))

    if a_ub is not None and b_ub is not None:
        ub_constraint_distance_vector = (b_ub - eps - a_ub @ current_location_tensor) / (a_ub @ unit_direction + eps)
        final_range_tensor = tf.minimum(final_range_tensor, filter_and_min(ub_constraint_distance_vector))
    elif a_ub is not None or b_ub is not None:
        raise ValueError("a_ub and b_ub should be both None or neither None!")
    success_tensor = final_range_tensor > eps
    final_range_tensor = tf.where(success_tensor, final_range_tensor, 0)
    return success_tensor, final_range_tensor


@tf_function(tf_function_mode)
def check_new_location(batch_size, new_location, lb, ub, a_ub, b_ub, eps):
    final_success_tensor = tf.constant([True] * batch_size)
    lb_check_tensor = tf.math.reduce_all(new_location - lb - eps >= 0, axis=1)
    final_success_tensor = tf.logical_and(final_success_tensor, lb_check_tensor)
    ub_check_tensor = tf.math.reduce_all(ub - new_location - eps >= 0, axis=1)
    final_success_tensor = tf.logical_and(final_success_tensor, ub_check_tensor)
    if a_ub is not None and b_ub is not None:
        ub_constraint_check_tensor = tf.math.reduce_all(b_ub - a_ub @ new_location - eps >= 0, axis=1)
        final_success_tensor = tf.logical_and(final_success_tensor, ub_constraint_check_tensor)
    return final_success_tensor


@tf_function(tf_function_mode)
def calculate_and_check_new_location(
        sampler_obj, current_location_tensor, direction_tensor, distance_1d_tensor, range_success_tensor, 
        batch_size):
    raw_next_location_tensor = current_location_tensor + direction_tensor * tf.reshape(
        distance_1d_tensor, [-1, 1])
    combined_next_location_tensor = tf.where(
        tf.reshape(range_success_tensor, [-1, 1]), raw_next_location_tensor, current_location_tensor)
    check_success_tensor = check_new_location(
        batch_size, combined_next_location_tensor, sampler_obj.min_value_vector, sampler_obj.max_value_vector,
        sampler_obj.smaller_eq_matrix, sampler_obj.smaller_eq_right_side_vector, sampler_obj.numeric_eps)
    final_success_tensor = tf.logical_and(range_success_tensor, check_success_tensor)
    return combined_next_location_tensor, final_success_tensor


@tf_function(tf_function_mode)
def update_location_and_check(
        self, batch_size, direction_range_tensor, current_location_tensor, combined_random_direction_tensor,
        range_success_tensor):
    random_distance = tensor_obj.random_uniform((batch_size,)) * direction_range_tensor * self.maximal_range_ratio
    combined_next_location_tensor, final_success_tensor = calculate_and_check_new_location(
        self, current_location_tensor, combined_random_direction_tensor, random_distance,
        range_success_tensor, batch_size)
    return combined_next_location_tensor, final_success_tensor


@tf_function(tf_function_mode)
def apply_location_update_and_count(
        self, final_success_tensor, combined_next_location_tensor, current_location_tensor,
        step_count_tensor, failed_count_tensor):
    final_next_location_tensor = tf.where(
        tf.reshape(final_success_tensor, [-1, 1]), combined_next_location_tensor, current_location_tensor)
    final_step_increment_tensor = tf.where(final_success_tensor, 1, 0)
    final_failed_assign_tensor = tf.where(final_success_tensor, 0, failed_count_tensor + 1)

    current_location_tensor.assign(final_next_location_tensor)
    step_count_tensor.assign_add(final_step_increment_tensor)
    failed_count_tensor.assign(final_failed_assign_tensor)

    updated_size = tf.math.count_nonzero(final_success_tensor, dtype=int_type)
    if updated_size > 0:
        updated_location_tensor = tf.boolean_mask(combined_next_location_tensor, final_success_tensor)
        self._update_center(updated_location_tensor, updated_size)
        self._update_saved_tensor(updated_location_tensor, updated_size)


@tf_function(tf_function_mode)
def drop_failed_and_add_new_locations(
        sampler_obj, failed_status_tensor, current_location_tensor, failed_count_tensor, step_count_tensor,
        batch_index_tensor=None, updated_counter_tensor=None):
    failed_num = tf.math.count_nonzero(failed_status_tensor, dtype=int_type)
    if failed_num > 0:
        failed_locations = tf.reshape(tf.where(failed_status_tensor), [-1])
        new_random_warmup_points = sampler_obj.random_warmup_point(failed_num)
        current_location_tensor.scatter_nd_update(failed_locations, new_random_warmup_points)
        final_failed_assign_tensor = tf.where(failed_status_tensor, 0, failed_count_tensor)
        failed_count_tensor.assign(final_failed_assign_tensor)
        final_step_count_assign_tensor = tf.where(failed_status_tensor, 0, step_count_tensor)
        step_count_tensor.assign(final_step_count_assign_tensor)
        if batch_index_tensor is not None:
            new_batch_index_tensor = tf.range(updated_counter_tensor, updated_counter_tensor + failed_num)
            updated_counter_tensor.assign_add(failed_num)
            batch_index_tensor.scatter_nd_update(failed_locations, new_batch_index_tensor)


def pickout_locations(pickout_status_tensor, current_location_tensor, pickout_tensor_storage_func):
    pickout_num = tf.math.count_nonzero(pickout_status_tensor, dtype=int_type)
    if pickout_num > 0:
        pickout_tensor = tf.boolean_mask(current_location_tensor, pickout_status_tensor)
        pickout_tensor_storage_func(pickout_tensor, pickout_num)


@tf_function(tf_function_mode)
def deal_with_pick_and_drop(
        self, failed_count_tensor, maximal_fail_num, step_count_tensor, maximal_sample_index,
        current_location_tensor, thinning_steps, pickout_tensor_storage_func):
    failed_status_tensor = tf.logical_or(
        failed_count_tensor > maximal_fail_num, step_count_tensor > maximal_sample_index)
    drop_failed_and_add_new_locations(
        self, failed_status_tensor, current_location_tensor, failed_count_tensor, step_count_tensor)

    pickout_status_tensor = tf.logical_and(step_count_tensor > 0, step_count_tensor % thinning_steps == 0)
    pickout_locations(pickout_status_tensor, current_location_tensor, pickout_tensor_storage_func)


@tf_function(tf_function_mode)
def new_step_function(
        self, current_location_tensor, step_count_tensor, failed_count_tensor, maximal_warmup_steps,
        maximal_sample_index, maximal_fail_num, thinning_steps, batch_size, pickout_tensor_storage_func):

    combined_random_direction_tensor = generate_random_direction(self, step_count_tensor, maximal_warmup_steps)

    range_success_tensor, direction_range_tensor = calculate_direction_range(
        batch_size, current_location_tensor, combined_random_direction_tensor, self.min_value_vector,
        self.max_value_vector, self.smaller_eq_matrix, self.smaller_eq_right_side_vector, self.numeric_eps)

    combined_next_location_tensor, final_success_tensor = update_location_and_check(
        self, batch_size, direction_range_tensor, current_location_tensor, combined_random_direction_tensor,
        range_success_tensor)

    apply_location_update_and_count(
        self, final_success_tensor, combined_next_location_tensor, current_location_tensor,
        step_count_tensor, failed_count_tensor)

    deal_with_pick_and_drop(
        self, failed_count_tensor, maximal_fail_num, step_count_tensor, maximal_sample_index,
        current_location_tensor, thinning_steps, pickout_tensor_storage_func)


@tf_function(tf_function_mode)
def sample_function(self, target_sample_size, output_sample_tensor, current_size):
    def pickout_tensor_storage_func(pickout_tensor, pickout_num):
        rest_num = target_sample_size - current_size
        if rest_num <= pickout_num:
            if rest_num < pickout_num:
                pickout_tensor = pickout_tensor[:rest_num]
                pickout_num = rest_num
        output_sample_tensor[current_size:current_size + pickout_num].assign(pickout_tensor)
        current_size.assign_add(pickout_num)

    each_step_iter_func = self.new_step
    while current_size < target_sample_size:
        each_step_iter_func(pickout_tensor_storage_func)

