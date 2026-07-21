from scripts.src.core.model.model_class import UserDefinedModel
from scripts.src.core.model.model_constructor import single_model_constructor, model_combine_constructor
from scripts.src.core.common.config import ModelKeyword

from .model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict


def tissue_specific_dict_generator(
        sink_tissue_list, source_tissue_list, base_model, serum_model, source_specific_model):
    tissue_specific_user_defined_model_dict = {
        ModelKeyword.serum: single_model_constructor(serum_model, model_metabolite_to_standard_name_dict)
    }
    base_reaction_list = []
    for pathway_name, pathway_reaction_list in base_model.reaction_dict.items():
        base_reaction_list.extend(pathway_reaction_list)
    for sink_tissue in sink_tissue_list:
        sink_user_defined_model = single_model_constructor(base_model, model_metabolite_to_standard_name_dict)
        tissue_specific_user_defined_model_dict[sink_tissue] = sink_user_defined_model
    for source_tissue in source_tissue_list:
        source_user_defined_model = model_combine_constructor(base_model, source_specific_model)
        tissue_specific_user_defined_model_dict[source_tissue] = source_user_defined_model
    return tissue_specific_user_defined_model_dict
