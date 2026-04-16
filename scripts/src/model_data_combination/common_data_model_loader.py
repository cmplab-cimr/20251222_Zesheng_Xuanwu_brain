from ..inventory import ExperimentName


def common_data_model_function_loader(model_name):
    if model_name == ExperimentName.data_20251115_ad_sema:
        from . import data_20251115_ad_sema as data_model_object
    elif model_name == ExperimentName.data_20251115_ad_sema_serum_labeled_buffer:
        from . import data_20251115_ad_sema_serum_labeled_buffer as data_model_object
    elif model_name == ExperimentName.data_20251115_ad_sema_single_tissue_model:
        from . import data_20251115_ad_sema_single_tissue_model as data_model_object
    else:
        raise ValueError()
    return data_model_object
