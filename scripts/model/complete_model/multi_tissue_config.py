from . import (invivo_infusion_model_with_tca_buffer, )

from scripts.src.core.common.config import ModelKeyword
from ..common_model_functions import single_model_constructor, model_combine_constructor
from ..model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict as standard_name_dict

complete_user_defined_infusion_model = single_model_constructor(invivo_infusion_model_with_tca_buffer, standard_name_dict)

gut_infusion_model_dict = {
    'gut': complete_user_defined_infusion_model,
}

whole_body_infusion_model_dict = {
    'whole_body': complete_user_defined_infusion_model,
}