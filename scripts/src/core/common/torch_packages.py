class TensorOps(object):
    def __init__(self, _float_type, _device):
        self.float_type = _float_type
        self.device = _device

    def ones(self, size, dtype=None, *args, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return torch.ones(size, dtype=dtype, device=self.device, *args, **kwargs)

    def zeros(self, size, dtype=None, *args, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return torch.zeros(size, dtype=dtype, device=self.device, *args, **kwargs)

    def zeros_like(self, target_tensor, *args, **kwargs):
        return torch.zeros_like(target_tensor, device=self.device, *args, **kwargs)

    def tensor(self, content, dtype=None, *args, **kwargs):
        if dtype is None:
            dtype = self.float_type
        return torch.tensor(content, dtype=dtype, device=self.device, *args, **kwargs)


try:
    import torch
except ModuleNotFoundError:
    torch = None
    torch_nn_func = None
    float_type = None
    device = None
    np_float_type = None
    int_type = None
    tensor_obj = None
else:
    from torch.nn import functional as torch_nn_func

    float_type = torch.double
    device = torch.device('cpu')
    # device = torch.device('cuda')
    np_float_type = 'float64'
    int_type = torch.long

    tensor_obj = TensorOps(float_type, device)

