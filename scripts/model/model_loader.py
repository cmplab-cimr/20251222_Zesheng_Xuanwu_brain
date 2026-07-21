from scripts.src.core.model.model_class import UserDefinedModel
from scripts.src.core.model.model_constructor import multiple_tissue_model_converter

from .model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict


class ModelList(object):
    liver_cortex_compartmentalized_model = 'liver_cortex_compartmentalized_model'
    liver_cortex_compartmentalized_model_with_serum_supplement = (
        'liver_cortex_compartmentalized_model_with_serum_supplement')
    liver_neuron_astrocyte_cortex_model = 'liver_neuron_astrocyte_cortex_model'
    tissue_specific_model = 'tissue_specific_model'


def multi_tissue_model_constructor(tissue_specific_model_dict):
    multiple_tissue_model = multiple_tissue_model_converter(tissue_specific_model_dict)
    return multiple_tissue_model


def model_loader(model_type, model_parameter=None):
    if model_type == ModelList.liver_cortex_compartmentalized_model:
        from .multi_tissue_compartmentalized_network.multi_tissue_config import liver_cortex_model_dict
        current_config = liver_cortex_model_dict
        current_model = multiple_tissue_model_converter(current_config)
    elif model_type == ModelList.liver_cortex_compartmentalized_model_with_serum_supplement:
        from .multi_tissue_compartmentalized_network.multi_tissue_config import liver_cortex_model_with_serum_labeled_buffer_generator
        current_config, current_labeled_buffer_input_metabolite_obj_dict = liver_cortex_model_with_serum_labeled_buffer_generator()
        current_user_defined_model = multiple_tissue_model_converter(current_config)
        current_model = (current_user_defined_model, current_labeled_buffer_input_metabolite_obj_dict)
    elif model_type == ModelList.tissue_specific_model:
        from .multi_tissue_compartmentalized_network.single_tissue_config import single_tissue_model_generator
        current_config, current_labeled_buffer_input_metabolite_obj_dict = single_tissue_model_generator(model_parameter)
        current_user_defined_model = multiple_tissue_model_converter(current_config)
        current_model = (current_user_defined_model, current_labeled_buffer_input_metabolite_obj_dict)
    elif model_type == ModelList.liver_neuron_astrocyte_cortex_model:
        from .multi_tissue_compartmentalized_network.multi_tissue_config import liver_complex_cortex_model_with_serum_buffer_generator
        current_config, *other_item = liver_complex_cortex_model_with_serum_buffer_generator()
        current_user_defined_model = multiple_tissue_model_converter(current_config)
        current_model = (current_user_defined_model, *other_item)
    else:
        raise ValueError()
    return current_model
