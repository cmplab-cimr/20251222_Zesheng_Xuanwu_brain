from scripts.src.core.model.model_class import UserDefinedModel
from scripts.src.core.model.model_constructor import multiple_tissue_model_converter

from .model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict


class ModelList(object):
    liver_cortex_compartmentalized_model = 'liver_cortex_compartmentalized_model'


def multi_tissue_model_constructor(tissue_specific_model_dict):
    multiple_tissue_model = multiple_tissue_model_converter(tissue_specific_model_dict)
    return multiple_tissue_model


def model_loader(model_type, model_parameter=None):
    if model_type == ModelList.liver_cortex_compartmentalized_model:
        from .multi_tissue_compartmentalized_network.multi_tissue_config import liver_cortex_model_dict
        current_config = liver_cortex_model_dict
    else:
        raise ValueError()
    current_model = multiple_tissue_model_converter(current_config)
    return current_model
