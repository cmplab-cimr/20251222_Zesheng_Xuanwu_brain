import os


class TensorOps(object):
    def __init__(self, _float_type):
        self.float_type = _float_type

    def ones(self, size, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.ones(size, *args, dtype=dtype, **kwargs)

    def ones_like(self, size, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.ones_like(size, *args, dtype=dtype, **kwargs)

    def zeros(self, size, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.zeros(size, *args, dtype=dtype, **kwargs)

    def zeros_like(self, target_tensor, *args, **kwargs):
        return tf.zeros_like(target_tensor, *args, **kwargs)

    def constant(self, content, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.constant(content, *args, dtype=dtype, **kwargs)

    def variable(self, initial_np_array, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.Variable(initial_np_array, *args, dtype=dtype, **kwargs)

    def variable_no_train(self, initial_np_array, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.Variable(initial_np_array, *args, dtype=dtype, trainable=False, **kwargs)

    def zero_variable_no_train(self, shape, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.Variable(zeros(shape, dtype=dtype), *args, dtype=dtype, trainable=False, **kwargs)

    def random_uniform(self, shape, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.random.uniform(shape, *args, dtype=dtype, **kwargs)

    def eye(self, size, *args, dtype=None, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return tf.eye(size, *args, dtype=dtype, **kwargs)

    def empty_sparse(self, complete_sparse_dim_shape, value_inner_shape, *args, **kwargs):
        non_zero_num = 0
        indices_tensor = self.zeros([non_zero_num], dtype=int_type)
        value_tensor = self.zeros([non_zero_num, *value_inner_shape], *args, **kwargs)
        new_empty_sparse = SparseTensor(
            indices_tensor, value_tensor, non_zero_num, complete_sparse_dim_shape, value_inner_shape)
        return new_empty_sparse


class SparseTensor(object):
    def __init__(self, indices_tensor, value_tensor, non_zero_num, sparse_dim_shape, value_inner_shape):
        self.indices_tensor = indices_tensor
        self.value_tensor = value_tensor
        self.non_zero_num = non_zero_num
        sparse_dim_shape = reshape(sparse_dim_shape, (-1,))
        self.sparse_dim_shape = sparse_dim_shape
        self.value_inner_shape = value_inner_shape
        self.complete_shape = tensor_obj.constant((*sparse_dim_shape, value_inner_shape), dtype=int_type)
        self.sparse_shape_cumulative_product = expand_dims(
            cumprod(self.sparse_dim_shape, exclusive=True, reverse=True), axis=-1)

    @staticmethod
    def construct_like(indices_tensor, value_tensor, non_zero_num, same_shape_sparse_tensor):
        return SparseTensor(
            indices_tensor, value_tensor, non_zero_num, same_shape_sparse_tensor.sparse_dim_shape,
            same_shape_sparse_tensor.value_inner_shape)

    def list_property(self):
        return (
            self.indices_tensor, self.value_tensor, self.non_zero_num, self.sparse_dim_shape,
            self.value_inner_shape, self.complete_shape)

    def flat_to_layer_index(self, flat_indices_tensor):
        return unravel_index(flat_indices_tensor, self.sparse_dim_shape)

    def layer_to_flat_index(self, layer_indices_tensor):
        if len(self.sparse_dim_shape) == 1:
            return layer_indices_tensor
        else:
            flat_index = reduce_sum(layer_indices_tensor * self.sparse_shape_cumulative_product, axis=1)
            return flat_index


# tf_function_mode = True
tf_function_mode = False
report_each_step = False
report_each_step = True
single_emu_calculation_test = False
efficient_float = True
xla_compile = None
if 'SINGLE_EMU_CALC_TEST' in os.environ:
    single_emu_calculation_test = True
    efficient_float = False
    # xla_compile = True
    xla_compile = None
# efficient_float = True


TensorSpec = None
tensor_array = None
ragged = None
transpose = None
device = None
float_type = None
int_type = None
uint_type = None
uint16 = None
uint8 = None
bool_type = None
tf_string_type = None
tensor_obj = None
constant = None
int_constant = None
zeros = None
variable = None

tf_function = None
tf_numpy_function = None
concat = None
qr_decomp = None
svd_decomp = None
eigh_decomp = None
tensor_pad = None
reshape = None
expand_dims = None
logical_and = None
argsort = None
argmax = None
sqrt = None
cholesky_decomp = None
cholesky_solve = None
triangular_solve = None
diag = None
diag_part = None
reverse = None
broadcast_to = None
floor = None
ceil = None
clip = None
type_conversion = None
minimum = None
maximum = None
boolean_mask = None
reduce_all = None
reduce_any = None
reduce_max = None
reduce_min = None
reduce_sum = None
unsorted_segment_sum = None
cumsum = None
scatter_nd = None
gather = None
gather_nd = None
stack = None
unstack = None
tensor_abs = None
is_nan = None
linear_solve = None
pinv = None
tensor_scatter_nd_update = None
tensor_scatter_nd_add = None
repeat = None
tile = None
reduce_mean = None
tensor_range = None
sign = None
ones_float = None
error_print = None
tf_print = None
transpose3d = None
count_nonzero = None
assign_and = None
assign_or = None
where = None
where_single = None
matrix_norm = None
least_square = None
argmin = None
cumprod = None
map_fn = None
lstsq = None
ragged_boolean_mask = None
tf_norm = None
tf_timestamp = None

inf = None

sparse_mat_mul = None

true_tensor = None
false_tensor = None
zero_1d_tensor = None
zero_2d_tensor = None
zero_3d_tensor = None
zero_4d_tensor = None
zero_5d_tensor = None
one_1d_tensor = None
one_2d_tensor = None
one_3d_tensor = None
inf_1d_tensor = None

ragged_constant = None
ragged_int_constant = None
ragged_range = None

report_step = None

try:
    import tensorflow as tf
except ModuleNotFoundError:
    tf = None
    tf_optimizers = None
else:
    from tensorflow import optimizers as tf_optimizers
    import functools
    import sys
    import numpy as np
    from sortedcontainers import SortedSet

    device = tf.device('CPU:0')
    # device = tf.device('GPU:0')
    float_type = tf.float64
    precise_float_type = tf.float64
    single_precision_float_type = tf.float32
    tf_string_type = tf.string
    if efficient_float:
        efficient_float_type = single_precision_float_type
    else:
        efficient_float_type = precise_float_type
    int_type = tf.int32
    int16 = tf.int16
    uint_type = tf.uint32
    uint16 = tf.uint16
    uint8 = tf.uint8
    bool_type = tf.bool
    TensorSpec = tf.TensorSpec

    tensor_obj = TensorOps(float_type)
    constant = tensor_obj.constant
    zeros = tensor_obj.zeros
    ones = tensor_obj.ones
    variable = tf.Variable

    def tf_function(mode=True, reduce_retracing=True, jit_compile=xla_compile, **kwargs):
        def identity(func):
            return func

        if mode:
            return tf.function(reduce_retracing=reduce_retracing, jit_compile=jit_compile, **kwargs)
        else:
            return identity

    def tf_numpy_function(*args, **kwargs):
        return tf.numpy_function(*args, **kwargs, Tout=float_type)

    def int_constant(*args, dtype=int_type, **kwargs):
        return tensor_obj.constant(*args, dtype=dtype, **kwargs)

    def tensor_array(*args, dtype=float_type, clear_after_read=True, **kwargs):
        return tf.TensorArray(dtype, size=0, dynamic_size=True, clear_after_read=clear_after_read)

    class TensorArray(object):
        def __init__(self, dtype=float_type, clear_after_read=True, **kwargs):
            self.array = tf.TensorArray(dtype, size=0, dynamic_size=True, clear_after_read=clear_after_read, **kwargs)
            self.count_index = 0

        def append(self, item_tensor):
            self.array.write(self.count_index, item_tensor)
            self.count_index += 1

    ragged = tf.ragged
    ragged_range = tf.ragged.range
    transpose = tf.transpose
    concat = tf.concat
    qr_decomp = tf.linalg.qr
    svd_decomp = tf.linalg.svd
    eigh_decomp = tf.linalg.eigh
    tensor_pad = tf.pad
    reshape = tf.reshape
    expand_dims = tf.expand_dims
    squeeze = tf.squeeze
    logical_and = tf.logical_and
    argmax = tf.math.argmax
    argsort = tf.argsort
    sqrt = tf.sqrt
    cholesky_decomp = tf.linalg.cholesky
    cholesky_solve = tf.linalg.cholesky_solve
    triangular_solve = tf.linalg.triangular_solve
    diag = tf.linalg.diag
    diag_part = tf.linalg.diag_part
    reverse = tf.reverse
    broadcast_to = tf.broadcast_to
    floor = tf.math.floor
    ceil = tf.math.ceil
    clip = tf.clip_by_value
    type_conversion = tf.cast
    minimum = tf.minimum
    maximum = tf.maximum
    boolean_mask = tf.boolean_mask
    reduce_all = tf.math.reduce_all
    reduce_any = tf.math.reduce_any
    reduce_max = tf.math.reduce_max
    reduce_min = tf.math.reduce_min
    reduce_sum = tf.math.reduce_sum
    reduce_mean = tf.reduce_mean
    unsorted_segment_sum = tf.math.unsorted_segment_sum
    cumsum = tf.math.cumsum
    scatter_nd = tf.scatter_nd
    gather = tf.gather
    gather_nd = tf.gather_nd
    stack = tf.stack
    unstack = tf.unstack
    tensor_abs = tf.abs
    is_nan = tf.math.is_nan
    linear_solve = tf.linalg.solve
    pinv = tf.linalg.pinv
    tensor_scatter_nd_update = tf.tensor_scatter_nd_update
    tensor_scatter_nd_add = tf.tensor_scatter_nd_add
    repeat = tf.repeat
    tile = tf.tile
    tensor_range = tf.range
    sign = tf.sign
    search_sorted = tf.searchsorted
    unravel_index = tf.unravel_index
    cumprod = tf.math.cumprod
    map_fn = tf.map_fn
    tf_print = tf.print
    tf_timestamp = tf.timestamp
    tf_str_format = tf.strings.format
    where = tf.where
    count_nonzero = tf.math.count_nonzero
    lstsq = tf.linalg.lstsq
    ragged_boolean_mask = tf.ragged.boolean_mask
    tf_norm = tf.norm

    sparse_mat_mul = tf.sparse.sparse_dense_matmul

    def ones_float(*args, **kwargs):
        return tf.ones(*args, dtype=float_type, **kwargs)

    @tf_function(True)
    def error_print(*args, **kwargs):
        tf.print(*args, output_stream=sys.stderr, **kwargs)

    def transpose3d(target_tensor, **kwargs):
        return tf.transpose(target_tensor, perm=[0, 2, 1], **kwargs)

    def assign_and_batch(variable, bool_tensor, b):
        variable[:b].assign(variable[:b] & bool_tensor)

    def assign_and(variable, bool_tensor):
        variable.assign(variable & bool_tensor)

    def assign_or(variable, bool_tensor):
        variable.assign(variable | bool_tensor)


    @tf_function(True)
    def where_single(condition):
        print('##### NNLS: Tracing where_single...')
        return type_conversion(tf.where(condition), dtype=int_type)

    def matrix_norm(a, *args, axis=None, **kwargs):
        if axis is None:
            axis = [-2, -1]
        return tf.norm(a, *args, axis=axis, **kwargs)

    def least_square(matrix, rhs):
        return tf.linalg.lstsq(matrix, rhs, fast=False)

    def argmin(*args, **kwargs):
        return tf.math.argmin(*args, output_type=int_type, **kwargs)

    inf = np.inf

    true_tensor = constant(True, dtype=bool_type)
    true_2d_tensor = constant([[True]], dtype=bool_type)
    false_tensor = constant(False, dtype=bool_type)
    zero_scalar_tensor = constant(0, dtype=precise_float_type)
    zero_1d_tensor = zeros(1)
    zero_2d_tensor = zeros((1, 1))
    zero_3d_tensor = zeros((1, 1, 1))
    zero_4d_tensor = zeros((1, 1, 1, 1))
    zero_5d_tensor = zeros((1, 1, 1, 1, 1))
    zero_scalar_int_tensor = constant(0, dtype=int_type)
    zero_1d_int_tensor = zeros(1, dtype=int_type)
    zero_2d_int_tensor = zeros((1, 1), dtype=int_type)
    one_1d_int_tensor = ones((1,), dtype=int_type)
    one_scalar_tensor = constant(1, dtype=precise_float_type)
    one_1d_tensor = ones(1)
    one_2d_tensor = ones((1, 1))
    one_3d_tensor = ones((1, 1, 1))
    inf_1d_tensor = ones(1) * inf
    one_3d_efficient_tensor = ones((1, 1, 1), dtype=efficient_float_type)
    one_1d_efficient_tensor = ones(1, dtype=efficient_float_type)
    zero_1d_efficient_tensor = zeros(1, dtype=efficient_float_type)
    zero_2d_efficient_tensor = zeros((1, 1), dtype=efficient_float_type)
    zero_3d_efficient_tensor = zeros((1, 1, 1), dtype=efficient_float_type)
    empty_str = constant('', dtype=tf_string_type)
    random_seed = tf.random.set_seed(4576251)
    shuffle = tf.random.shuffle
    tf_uniform = tf.random.uniform

    special_report_output_stream = sys.stderr
    normal_report_output_stream = sys.stdout

    perm_3d = int_constant((0, 2, 1))

    if report_each_step:
        report_step = true_tensor
    else:
        report_step = false_tensor

    def ragged_constant(*args, dtype=float_type, **kwargs):
        return ragged.constant(*args, dtype=dtype, **kwargs)

    def ragged_int_constant(*args, dtype=int_type, **kwargs):
        return ragged.constant(*args, dtype=dtype, **kwargs)

    @tf_function(True)
    def fake_timestamp():
        print('!!!!: Tracing fake_timestamp...')
        return zero_1d_tensor

    @tf_function(True)
    def fake_str_format(template, inputs):
        return constant('', dtype=tf_string_type)

    @tf_function(True)
    def format_time(float_time):
        print('!!!!: Tracing format_time...')
        hour_offset = -4
        time_stamp = type_conversion(float_time, int_type)
        hour_minute_second = time_stamp % 86400
        hour = (hour_minute_second // 3600 + hour_offset) % 24
        minute_second = hour_minute_second % 3600
        minute = minute_second // 60
        second = minute_second % 60
        return hour, minute, second

    @tf_function(True)
    def format_time_with_millisecond(float_time):
        print('!!!!: Tracing format_time_with_millisecond...')
        hour, minute, second = format_time(float_time)
        millisecond = type_conversion(float_time % 1 * 1000, int_type) % 1000
        return hour, minute, second, millisecond

    @tf_function(True)
    def normal_output_with_time(formatted_str, float_time, special_output):
        print('!!!!: Tracing normal_output_with_time...')
        hour, minute, second, millisecond = format_time_with_millisecond(float_time)
        time_str = tf.strings.format('[{}:{}:{}.{}]', (hour, minute, second, millisecond))
        final_str = tf.strings.join([time_str, formatted_str], separator=' ')
        if special_output:
            tf_print(final_str, output_stream=sys.stderr)
        else:
            tf_print(final_str, output_stream=sys.stdout)
        return None

    @tf_function(True)
    def fake_output_with_time(formatted_str, float_time, special_output):
        print('!!!!: Tracing fake_output_with_time...')
        return None

    @tf_function(True)
    def fake_format_and_output_with_time(template_index, inputs, float_time, special_output):
        print('!!!!: Tracing fake_format_and_output_with_time...')
        return None

    def test_map_fn(func, parameter_tuple, parallel_iterations, fn_output_signature, n=0):
        if isinstance(parameter_tuple, tuple):
            input_num = len(parameter_tuple)
            batch_size = parameter_tuple[0].shape[0]
        else:
            input_num = 1
            batch_size = parameter_tuple.shape[0]
        result_num = len(fn_output_signature)
        output_list = [[] for _ in range(result_num)]
        for i in range(batch_size):
            if input_num == 1:
                arg_list = parameter_tuple[i]
            else:
                arg_list = [parameter_tuple[j][i] for j in range(input_num)]
            returned_value_list = func(arg_list)
            if result_num == 1:
                output_list[0].append(returned_value_list)
            else:
                for k in range(result_num):
                    output_list[k].append(returned_value_list[k])
        output_tensor_list = [stack(output_item) for output_item in output_list]
        return output_tensor_list

    @tf_function(True)
    def while_map_func(fn, elems, parallel_iterations, fn_output_signature, n):
        def compute(i, result_tensor_array_list):
            print('############### Tracing while_map compute...')
            result_tuple = fn([each_element[i] for each_element in elems])
            result_tensor_array_list = [
                ta.write(i, value) for (ta, value) in zip(result_tensor_array_list, result_tuple)
            ]
            return i + 1, result_tensor_array_list

        result_tensor_array_list = [
            tf.TensorArray(dtype=output_dtype, size=n, dynamic_size=False)
            for output_dtype in fn_output_signature]
        for i in range(n):
            result_tuple = fn([each_element[i] for each_element in elems])
            result_tensor_array_list = [
                ta.write(i, value) for (ta, value) in zip(result_tensor_array_list, result_tuple)]

        # final_tensor_array_list = tf.while_loop(
        #     lambda i, _: i < n,
        #     compute, (initial_i, result_tensor_array_list),
        #     parallel_iterations=parallel_iterations,
        #     maximum_iterations=n)
        return [final_tensor_array.stack() for final_tensor_array in result_tensor_array_list]


def batch_no_train_variable_maker(shape_or_initial, dtype=None, total_num=1, zero=False, **kwargs):
    if zero:
        zero_initial_tensor = zeros(shape_or_initial, dtype=dtype)
        return [tensor_obj.variable_no_train(zero_initial_tensor, dtype=dtype, **kwargs) for _ in range(total_num)]
    else:
        return [tensor_obj.variable_no_train(shape_or_initial, dtype=dtype, **kwargs) for _ in range(total_num)]

