from scripts.src.core.model.model_class import UserDefinedModel
from scripts.src.core.model.model_constructor import multiple_tissue_model_converter
from scripts.src.core.common.config import ModelKeyword

from .model_metabolite_to_standard_name_dict import model_metabolite_to_standard_name_dict


def check_attribute(model_obj):
    try:
        model_obj.symmetrical_metabolite_set
    except AttributeError:
        model_obj.symmetrical_metabolite_set = set()
    try:
        model_obj.added_input_metabolite_set
    except AttributeError:
        model_obj.added_input_metabolite_set = set()
    try:
        model_obj.emu_excluded_metabolite_set
    except AttributeError:
        model_obj.emu_excluded_metabolite_set = set()
    try:
        model_obj.balance_excluded_metabolite_set
    except AttributeError:
        model_obj.balance_excluded_metabolite_set = set()
    try:
        model_obj.model_compartment_set
    except AttributeError:
        model_obj.model_compartment_set = set()
    try:
        model_obj.composite_reaction_list
    except AttributeError:
        model_obj.composite_reaction_list = []
    try:
        model_obj.other_pool_metabolite_dict
    except AttributeError:
        model_obj.other_pool_metabolite_dict = {
                ModelKeyword.general: model_obj.emu_excluded_metabolite_set}


def single_model_constructor(base_model, _model_metabolite_to_standard_name_dict):
    check_attribute(base_model)
    return UserDefinedModel(
            reaction_list=base_model.reaction_list,
            symmetrical_metabolite_set=base_model.symmetrical_metabolite_set,
            added_input_metabolite_set=base_model.added_input_metabolite_set,
            emu_excluded_metabolite_set=base_model.emu_excluded_metabolite_set,
            balance_excluded_metabolite_set=base_model.balance_excluded_metabolite_set,
            target_metabolite_list=None,
            model_compartment_set=base_model.model_compartment_set,
            composite_reaction_list=base_model.composite_reaction_list,
            other_pool_metabolite_dict=base_model.other_pool_metabolite_dict,
            model_metabolite_to_standard_name_dict=_model_metabolite_to_standard_name_dict,)


def model_combine_constructor(base_model, specific_addon_model, _model_metabolite_to_standard_name_dict):
    check_attribute(base_model)
    check_attribute(specific_addon_model)
    source_reaction_list = base_model.reaction_list + specific_addon_model.extra_reaction_list
    source_emu_excluded_metabolite_set = base_model.emu_excluded_metabolite_set | \
        specific_addon_model.emu_excluded_metabolite_set
    source_balance_excluded_metabolite_set = base_model.balance_excluded_metabolite_set | \
        specific_addon_model.balance_excluded_metabolite_set
    source_symmetrical_metabolite_set = base_model.symmetrical_metabolite_set | \
        specific_addon_model.symmetrical_metabolite_set
    source_added_input_metabolite_set = base_model.added_input_metabolite_set | \
        specific_addon_model.added_input_metabolite_set
    source_model_compartment_set = base_model.model_compartment_set | specific_addon_model.model_compartment_set
    source_composite_reaction_list = base_model.composite_reaction_list + \
        specific_addon_model.composite_reaction_list
    source_other_pool_metabolite_dict = {}
    current_other_pool_metabolite_dict = specific_addon_model.other_pool_metabolite_dict
    for tissue_key, raw_set in base_model.other_pool_metabolite_dict.items():
        new_set = set(raw_set)
        if tissue_key in current_other_pool_metabolite_dict:
            new_set |= current_other_pool_metabolite_dict[tissue_key]
        source_other_pool_metabolite_dict[tissue_key] = new_set
    source_user_defined_model = UserDefinedModel(
        reaction_list=source_reaction_list,
        symmetrical_metabolite_set=source_symmetrical_metabolite_set,
        added_input_metabolite_set=source_added_input_metabolite_set,
        emu_excluded_metabolite_set=source_emu_excluded_metabolite_set,
        balance_excluded_metabolite_set=source_balance_excluded_metabolite_set,
        target_metabolite_list=None,
        model_compartment_set=source_model_compartment_set,
        composite_reaction_list=source_composite_reaction_list,
        other_pool_metabolite_dict=source_other_pool_metabolite_dict,
        model_metabolite_to_standard_name_dict=_model_metabolite_to_standard_name_dict,
    )
    return source_user_defined_model


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
