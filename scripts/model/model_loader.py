from scripts.src.core.model.model_class import UserDefinedModel
from scripts.src.core.model.model_constructor import multiple_tissue_model_converter

from .model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict
from .complete_model.multi_tissue_config import (
    gut_infusion_model_dict, whole_body_infusion_model_dict)


class ModelList(object):
    gut_infusion_model_dict = 'gut_infusion_model_dict'
    whole_body_infusion_model_dict = 'whole_body_infusion_model_dict'


def multi_tissue_model_constructor(tissue_specific_model_dict):
    multiple_tissue_model = multiple_tissue_model_converter(tissue_specific_model_dict)
    return multiple_tissue_model


def model_loader(model_type, model_parameter=None):
    if model_type == ModelList.gut_infusion_model_dict:
        current_config = gut_infusion_model_dict
    elif model_type == ModelList.whole_body_infusion_model_dict:
        current_config = whole_body_infusion_model_dict
    else:
        raise ValueError()
    current_model = multiple_tissue_model_converter(current_config)
    return current_model
