from . import (
    basic_tissue_compartmentalized_model, liver_specific_compartmentalized_model,
    ketone_body_utilization_specific_compartmentalized_model, serum_specific_model, muscle_specific_model,
    adipose_specific_model)

from scripts.src.core.common.config import ModelKeyword
from ..common_model_functions import single_model_constructor, model_combine_constructor
from ..model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict as standard_name_dict

typical_liver_model = model_combine_constructor(
    basic_tissue_compartmentalized_model, liver_specific_compartmentalized_model, standard_name_dict)
typical_brain_model = model_combine_constructor(
    basic_tissue_compartmentalized_model, ketone_body_utilization_specific_compartmentalized_model, standard_name_dict)

liver_cortex_model_dict = {
    ModelKeyword.serum: single_model_constructor(serum_specific_model, standard_name_dict),

    ModelKeyword.liver: typical_liver_model,
    ModelKeyword.cortex: typical_brain_model,
}
