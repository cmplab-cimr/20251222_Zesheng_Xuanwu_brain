from ..inventory import ExperimentName


def common_data_model_function_loader(model_name):
    if model_name == ExperimentName.example_data_20250815:
        from . import example_data_20250815 as data_model_object
    else:
        raise ValueError()
    return data_model_object
