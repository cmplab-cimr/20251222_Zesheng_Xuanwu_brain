from ..common.config import ModelKeyword
from ..common.functions import tissue_specific_name_constructor, reverse_reaction_name
from ..common.classes import TransformDict

from .model_class import EMUMIDDimDict, MFAModel, UserDefinedModel
from .common_model_process_functions import model_preprocess, compart_all_metabolites, \
    flux_balance_equation_generator, metabolite_carbon_number_verification
from .emu_analyzer_functions import emu_equation_analyzer, emu_matrix_equation_generator, emu_dependency_analyzer


def common_model_constructor(user_defined_model: UserDefinedModel, new_format=False):
    (
        metabolite_reaction_dict, product_reaction_dict_for_emu, composite_reaction_dict,
        flux_name_index_dict, input_metabolite_name_set, complete_metabolite_dim_dict,
        target_metabolite_name_list) = model_preprocess(
        user_defined_model.reaction_list, user_defined_model.symmetrical_metabolite_set,
        user_defined_model.added_input_metabolite_set, user_defined_model.emu_excluded_metabolite_set,
        user_defined_model.balance_excluded_metabolite_set, user_defined_model.target_metabolite_list,
        user_defined_model.composite_reaction_list)

    all_target_metabolite_name_set = list(complete_metabolite_dim_dict.keys() - input_metabolite_name_set)
    # all_target_metabolite_name_set = {
    #     key: value for key, value in complete_metabolite_dim_dict.items() if key not in input_metabolite_name_set}
    (
        complete_tissue_metabolite_compartment_dict, metabolite_bare_metabolite_name_dict,
        metabolite_full_name_tissue_compartment_dict
    ) = compart_all_metabolites(
        complete_metabolite_dim_dict.keys(), user_defined_model.model_compartment_set)
    bare_metabolite_dim_dict = metabolite_carbon_number_verification(
        complete_metabolite_dim_dict, metabolite_bare_metabolite_name_dict)
    flux_balance_matrix, flux_balance_right_side_vector = flux_balance_equation_generator(
        metabolite_reaction_dict, composite_reaction_dict, flux_name_index_dict)

    emu_mid_equation_dict, input_emu_dict = emu_equation_analyzer(
        product_reaction_dict_for_emu, target_metabolite_name_list, input_metabolite_name_set,
        complete_metabolite_dim_dict)
    output_emu_mid_equation_dict = emu_matrix_equation_generator(emu_mid_equation_dict, input_emu_dict)
    emu_name_dependency_dict, complete_emu_obj_index_dict, new_input_emu_dict = emu_dependency_analyzer(
        product_reaction_dict_for_emu, input_metabolite_name_set, complete_metabolite_dim_dict)
    assert input_emu_dict.keys() == new_input_emu_dict.keys()
    complete_emu_dim_dict = EMUMIDDimDict()
    mfa_model = MFAModel(
        input_emu_dict, target_metabolite_name_list, output_emu_mid_equation_dict, emu_name_dependency_dict,
        composite_reaction_dict, complete_emu_dim_dict, complete_emu_obj_index_dict, flux_name_index_dict,
        flux_balance_matrix, flux_balance_right_side_vector, bare_metabolite_dim_dict,
        metabolite_bare_metabolite_name_dict, complete_tissue_metabolite_compartment_dict,
        input_metabolite_name_set, all_target_metabolite_name_set,
        user_defined_model.model_metabolite_to_standard_name_dict, user_defined_model,
        metabolite_full_name_tissue_compartment_dict)
    return mfa_model


def multiple_tissue_model_converter(tissue_specific_user_defined_model_dict):
    def modified_tissue_name_generator(raw_reactant_name, current_tissue_name, other_compartment_metabolite_dict):
        if raw_reactant_name in other_compartment_metabolite_dict:
            target_tissue_name = other_compartment_metabolite_dict[raw_reactant_name]
            if target_tissue_name == ModelKeyword.general:
                modified_tissue_name = None
            else:
                modified_tissue_name = target_tissue_name
        else:
            modified_tissue_name = current_tissue_name
        return modified_tissue_name

    def rename_reaction_list(raw_reaction_list, current_tissue_name, other_compartment_metabolite_dict):
        new_reaction_list = []
        for raw_reaction in raw_reaction_list:
            raw_reactant_name = raw_reaction[0]
            modified_tissue_name = modified_tissue_name_generator(
                raw_reactant_name, current_tissue_name, other_compartment_metabolite_dict)
            if modified_tissue_name is None:
                new_reaction_list.append(raw_reaction)
                continue
            new_reaction = list(raw_reaction)
            new_reaction[0] = tissue_specific_name_constructor(raw_reactant_name, modified_tissue_name)
            new_reaction_list.append(tuple(new_reaction))
        return new_reaction_list

    def rename_and_add_to_complete_set(
            raw_single_set, complete_set, current_tissue_name, other_compartment_metabolite_dict):
        for raw_name in raw_single_set:
            modified_tissue_name = modified_tissue_name_generator(
                raw_name, current_tissue_name, other_compartment_metabolite_dict)
            new_name = tissue_specific_name_constructor(raw_name, modified_tissue_name)
            complete_set.add(new_name)

    complete_reaction_list = []
    complete_emu_excluded_metabolite_set = set()
    complete_balance_excluded_metabolite_set = set()
    complete_symmetrical_metabolite_set = set()
    complete_added_input_metabolite_set = set()
    complete_model_compartment_set = set()
    complete_model_metabolite_to_standard_name_dict = TransformDict()
    complete_model_tissue_set = set(tissue_specific_user_defined_model_dict.keys())
    complete_composite_reaction_list = []
    for tissue_name, each_tissue_user_defined_model in tissue_specific_user_defined_model_dict.items():
        current_other_compartment_metabolite_dict = {}
        for other_tissue_name, other_metabolite_set in \
                each_tissue_user_defined_model.other_pool_metabolite_dict.items():
            for other_metabolite_name in other_metabolite_set:
                current_other_compartment_metabolite_dict[other_metabolite_name] = other_tissue_name
        for raw_reaction_dict in each_tissue_user_defined_model.reaction_list:
            new_reaction_dict = dict(raw_reaction_dict)
            raw_id = new_reaction_dict[ModelKeyword.id]
            raw_sub_reaction_list = new_reaction_dict[ModelKeyword.sub]
            raw_pro_reaction_list = new_reaction_dict[ModelKeyword.pro]
            new_id = tissue_specific_name_constructor(raw_id, tissue_name)
            new_sub_reaction_list = rename_reaction_list(
                raw_sub_reaction_list, tissue_name, current_other_compartment_metabolite_dict)
            new_pro_reaction_list = rename_reaction_list(
                raw_pro_reaction_list, tissue_name, current_other_compartment_metabolite_dict)
            new_reaction_dict.update({
                ModelKeyword.id: new_id,
                ModelKeyword.sub: new_sub_reaction_list,
                ModelKeyword.pro: new_pro_reaction_list})
            complete_reaction_list.append(new_reaction_dict)

        raw_composite_reaction_list = each_tissue_user_defined_model.composite_reaction_list
        if raw_composite_reaction_list is not None:
            for each_composite_reaction_dict in raw_composite_reaction_list:
                raw_id = each_composite_reaction_dict[ModelKeyword.id]
                raw_comp_list = each_composite_reaction_dict[ModelKeyword.comp]
                new_id = tissue_specific_name_constructor(raw_id, tissue_name)
                new_comp_list = rename_reaction_list(
                    raw_comp_list, tissue_name, current_other_compartment_metabolite_dict)
                new_composite_reaction_dict = {
                    **each_composite_reaction_dict,
                    ModelKeyword.id: new_id,
                    ModelKeyword.comp: new_comp_list,
                }
                complete_composite_reaction_list.append(new_composite_reaction_dict)

        raw_emu_excluded_metabolite_set = each_tissue_user_defined_model.emu_excluded_metabolite_set
        rename_and_add_to_complete_set(
            raw_emu_excluded_metabolite_set, complete_emu_excluded_metabolite_set, tissue_name,
            current_other_compartment_metabolite_dict)
        raw_balance_excluded_metabolite_set = each_tissue_user_defined_model.balance_excluded_metabolite_set
        rename_and_add_to_complete_set(
            raw_balance_excluded_metabolite_set, complete_balance_excluded_metabolite_set, tissue_name,
            current_other_compartment_metabolite_dict)
        raw_symmetrical_metabolite_set = each_tissue_user_defined_model.symmetrical_metabolite_set
        rename_and_add_to_complete_set(
            raw_symmetrical_metabolite_set, complete_symmetrical_metabolite_set, tissue_name,
            current_other_compartment_metabolite_dict)
        raw_added_input_metabolite_set = each_tissue_user_defined_model.added_input_metabolite_set
        rename_and_add_to_complete_set(
            raw_added_input_metabolite_set, complete_added_input_metabolite_set, tissue_name,
            current_other_compartment_metabolite_dict)
        raw_model_compartment_set = each_tissue_user_defined_model.model_compartment_set
        complete_model_compartment_set |= raw_model_compartment_set
        model_metabolite_to_standard_name_dict = each_tissue_user_defined_model.model_metabolite_to_standard_name_dict
        complete_model_metabolite_to_standard_name_dict.update(model_metabolite_to_standard_name_dict)

    complete_user_defined_model = UserDefinedModel(
        reaction_list=complete_reaction_list,
        symmetrical_metabolite_set=complete_symmetrical_metabolite_set,
        added_input_metabolite_set=complete_added_input_metabolite_set,
        emu_excluded_metabolite_set=complete_emu_excluded_metabolite_set,
        balance_excluded_metabolite_set=complete_balance_excluded_metabolite_set,
        target_metabolite_list=None,
        model_compartment_set=complete_model_compartment_set,
        model_tissue_set=complete_model_tissue_set,
        composite_reaction_list=complete_composite_reaction_list,
        model_metabolite_to_standard_name_dict=complete_model_metabolite_to_standard_name_dict,
    )
    return complete_user_defined_model


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


def remove_reaction(source_reaction_list, removed_reaction_id_pair_dict):
    forward_label = 'F'
    reverse_label = 'R'
    bidirectional_label = 'B'
    rxn_reversibility_label = 'reverse'
    rxn_product_label = 'pro'
    rxn_substrate_label = 'sub'
    rxn_id_label = 'id'

    present_ids = {rxn[rxn_id_label] for rxn in source_reaction_list}
    unmatched_id_set = set()
    for removed_id, remove_direction in removed_reaction_id_pair_dict.items():
        if removed_id not in present_ids:
            unmatched_id_set.add(removed_id)
    if len(unmatched_id_set) > 0:
        import warnings
        warnings.warn(
            f"model_combine_constructor: 下列待删反应 id 在合并模型中不存在,已忽略: {sorted(unmatched_id_set)}")
    new_source_reaction_list = []
    for rxn in source_reaction_list:
        rxn_id = rxn[rxn_id_label]
        if rxn_id in removed_reaction_id_pair_dict:
            remove_direction = removed_reaction_id_pair_dict[rxn_id]
            reversibility = rxn[rxn_reversibility_label] if rxn_reversibility_label in rxn else False
            new_rxn = None
            if not reversibility:
                assert remove_direction == forward_label
            else:
                if remove_direction == bidirectional_label:
                    pass
                elif remove_direction == reverse_label:
                    new_rxn = {
                        rxn_label: rxn_content for rxn_label, rxn_content in rxn.items()
                        if rxn_label != rxn_reversibility_label}
                elif remove_direction == forward_label:
                    new_rxn = {
                        k: v for k, v in rxn.items()
                        if k != rxn_reversibility_label}
                    new_rxn[rxn_id_label] = reverse_reaction_name(rxn_id)
                    new_rxn[rxn_product_label], new_rxn[rxn_substrate_label] = (
                        rxn[rxn_substrate_label], rxn[rxn_product_label])
                else:
                    raise ValueError()
            if new_rxn is not None:
                new_source_reaction_list.append(new_rxn)
        else:
            new_source_reaction_list.append(rxn)
    return new_source_reaction_list


def model_combine_constructor(
        base_model, specific_addon_model,
        _model_metabolite_to_standard_name_dict,
        removed_reaction_id_pair_dict = None,
        removed_emu_excluded_metabolite_set = None,
        removed_balance_excluded_metabolite_set = None,
        removed_symmetrical_metabolite_set = None,
        removed_added_input_metabolite_set = None,
        removed_model_compartment_set = None,
        removed_other_pool_metabolite_dict = None,
):
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
    for tissue_key, raw_set in current_other_pool_metabolite_dict.items():
        if tissue_key not in source_other_pool_metabolite_dict:
            source_other_pool_metabolite_dict[tissue_key] = set(raw_set)

    if removed_reaction_id_pair_dict:
        source_reaction_list = remove_reaction(source_reaction_list, removed_reaction_id_pair_dict)
    if removed_emu_excluded_metabolite_set:
        source_emu_excluded_metabolite_set = \
            source_emu_excluded_metabolite_set - set(removed_emu_excluded_metabolite_set)
    if removed_balance_excluded_metabolite_set:
        source_balance_excluded_metabolite_set = \
            source_balance_excluded_metabolite_set - set(removed_balance_excluded_metabolite_set)
    if removed_symmetrical_metabolite_set:
        source_symmetrical_metabolite_set = \
            source_symmetrical_metabolite_set - set(removed_symmetrical_metabolite_set)
    if removed_added_input_metabolite_set:
        source_added_input_metabolite_set = \
            source_added_input_metabolite_set - set(removed_added_input_metabolite_set)
    if removed_model_compartment_set:
        source_model_compartment_set = \
            source_model_compartment_set - set(removed_model_compartment_set)
    if removed_other_pool_metabolite_dict:
        for tissue_key, met_set in removed_other_pool_metabolite_dict.items():
            if tissue_key in source_other_pool_metabolite_dict:
                source_other_pool_metabolite_dict[tissue_key] = \
                    source_other_pool_metabolite_dict[tissue_key] - set(met_set)
            else:
                import warnings
                warnings.warn(
                    f"model_combine_constructor: 下列待删代谢物组织 tissue_key "
                    f"在合并模型中不存在,已忽略: {tissue_key}: {str(met_set)}")

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
