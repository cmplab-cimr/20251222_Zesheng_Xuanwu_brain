from . import invivo_infusion_model
from scripts.src.core.common.config import ModelKeyword

reaction_list = [
    *invivo_infusion_model.reaction_list,
    {
        'id': 'PYR_supplement',
        'sub': [('PYR_stock', 'abc')],
        'pro': [('PYR_c', 'abc')],
        'reverse': True
    },
    {
        'id': 'CIT_supplement',
        'sub': [('CIT_stock', 'abcdef')],
        'pro': [('CIT_m', 'abcdef')],
        'reverse': True
    },
]

emu_excluded_metabolite_set = {
    *invivo_infusion_model.emu_excluded_metabolite_set,
    'PYR_stock',
    'CIT_stock',
}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

symmetrical_metabolite_set = invivo_infusion_model.symmetrical_metabolite_set

added_input_metabolite_set = invivo_infusion_model.added_input_metabolite_set

model_compartment_set = invivo_infusion_model.model_compartment_set

composite_reaction_list = [
    *invivo_infusion_model.composite_reaction_list,
    {
        'id': 'PYR_supplement_net',
        'comp': [('PYR_supplement', ), ('PYR_supplement__R', -1)],
        ModelKeyword.flux_range: ModelKeyword.add_range_type
    },
    {
        'id': 'CIT_supplement_net',
        'comp': [('CIT_supplement', ), ('CIT_supplement__R', -1)],
        ModelKeyword.flux_range: ModelKeyword.add_range_type
    },
]

