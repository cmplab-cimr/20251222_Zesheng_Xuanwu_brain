from . import (
    basic_tissue_compartmentalized_model, liver_specific_compartmentalized_model,
    ketone_body_utilization_specific_compartmentalized_model, serum_specific_model, muscle_specific_model,
    adipose_specific_model, kidney_specific_compartmentalized_model)

from scripts.src.core.common.config import ModelKeyword
from ..common_model_functions import single_model_constructor, model_combine_constructor
from ..model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict as standard_name_dict


def single_tissue_model_generator(model_parameter_dict):
    from . import serum_metabolite_labeled_supplement_model
    from scripts.src.core.data.data_class import input_metabolite_obj_dict_batch_generator
    labeled_buffer_metabolite_name_carbon_num_list = (
        serum_metabolite_labeled_supplement_model.labeled_buffer_metabolite_name_carbon_num_list)
    labeled_buffer_input_metabolite_obj_dict = input_metabolite_obj_dict_batch_generator(
        labeled_buffer_metabolite_name_carbon_num_list)
    current_tissue = model_parameter_dict['tissue']
    serum_model = model_combine_constructor(
        serum_specific_model, serum_metabolite_labeled_supplement_model, standard_name_dict)
    complete_model_dict = {
        ModelKeyword.serum: serum_model,}
    if current_tissue == ModelKeyword.liver:
        tissue_model = model_combine_constructor(
            basic_tissue_compartmentalized_model,
            liver_specific_compartmentalized_model, standard_name_dict)
    elif current_tissue == ModelKeyword.kidney:
        tissue_model = model_combine_constructor(
            basic_tissue_compartmentalized_model,
            kidney_specific_compartmentalized_model, standard_name_dict)
    elif current_tissue in {
        ModelKeyword.cortex, ModelKeyword.hippocampus, ModelKeyword.hypothalamus}:
        tissue_model = model_combine_constructor(
            basic_tissue_compartmentalized_model,
            ketone_body_utilization_specific_compartmentalized_model, standard_name_dict)
    elif current_tissue in {
        ModelKeyword.q_muscle, ModelKeyword.soleus_muscle
    }:
        tissue_model = model_combine_constructor(
            model_combine_constructor(
                basic_tissue_compartmentalized_model,
                ketone_body_utilization_specific_compartmentalized_model,
                standard_name_dict),
            muscle_specific_model, standard_name_dict)
    elif current_tissue == ModelKeyword.brown_adipose:
        tissue_model = model_combine_constructor(
            model_combine_constructor(
                basic_tissue_compartmentalized_model, adipose_specific_model, standard_name_dict),
            ketone_body_utilization_specific_compartmentalized_model, standard_name_dict)
    elif current_tissue == ModelKeyword.colon:
        tissue_model = model_combine_constructor(
            basic_tissue_compartmentalized_model, ketone_body_utilization_specific_compartmentalized_model, standard_name_dict)
    else:
        raise ValueError()
    complete_model_dict.update({
        current_tissue: tissue_model
    })
    return complete_model_dict, labeled_buffer_input_metabolite_obj_dict
