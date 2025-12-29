from . import (
    basic_tissue_compartmentalized_model, liver_specific_compartmentalized_model,
    ketone_body_utilization_specific_compartmentalized_model, serum_specific_model, muscle_specific_model,
    adipose_specific_model)

from scripts.src.core.common.config import ModelKeyword
from ..common_model_functions import single_model_constructor, model_combine_constructor
from ..model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict as standard_name_dict


def serum_model_constructor(with_sink_fluxes=True):
    if with_sink_fluxes:
        from . import serum_metabolite_buffer_model
        serum_model = model_combine_constructor(
            serum_specific_model, serum_metabolite_buffer_model, standard_name_dict)
    else:
        serum_model = single_model_constructor(serum_specific_model, standard_name_dict)
    return serum_model


def typical_liver_model_constructor():
    typical_liver_model = model_combine_constructor(
        basic_tissue_compartmentalized_model, liver_specific_compartmentalized_model, standard_name_dict)
    return typical_liver_model


def typical_brain_model_constructor():
    typical_brain_model = model_combine_constructor(
        basic_tissue_compartmentalized_model, ketone_body_utilization_specific_compartmentalized_model, standard_name_dict)
    return typical_brain_model

liver_cortex_model_dict = {
    ModelKeyword.serum: serum_model_constructor(with_sink_fluxes=True),

    ModelKeyword.liver: typical_liver_model_constructor(),
    ModelKeyword.cortex: typical_brain_model_constructor(),
    # ModelKeyword.kidney: typical_liver_model_constructor(),
    # ModelKeyword.brown_adipose: model_combine_constructor(
    #     model_combine_constructor(
    #         basic_tissue_compartmentalized_model, adipose_specific_model, standard_name_dict),
    #     ketone_body_utilization_specific_compartmentalized_model, standard_name_dict),
}
