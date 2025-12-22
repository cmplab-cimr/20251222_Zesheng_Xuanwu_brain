from . import complete_model
from scripts.src.core.common.config import ModelKeyword

reaction_list = [
    *complete_model.reaction_list,
    {
        'id': 'GLC_unlabelled_input',
        'sub': [('GLC_unlabelled_e', 'abcdef')],
        'pro': [('GLC_c', 'abcdef')],
    },
    # {
    #     'id': 'PYR_unlabelled_input',
    #     'sub': [('PYR_unlabelled_c', 'abc')],
    #     'pro': [('PYR_c', 'abc')],
    # },
    # {
    #     'id': 'CIT_unlabelled_input',
    #     'sub': [('CIT_unlabelled_m', 'abcdef')],
    #     'pro': [('CIT_m', 'abcdef')],
    # },
]

emu_excluded_metabolite_set = {
    *complete_model.emu_excluded_metabolite_set,
    'GLC_unlabelled_e',
    # 'PYR_unlabelled_c',
    # 'CIT_unlabelled_m',
}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

symmetrical_metabolite_set = complete_model.symmetrical_metabolite_set

added_input_metabolite_set = complete_model.added_input_metabolite_set

model_compartment_set = complete_model.model_compartment_set

composite_reaction_list = [
    {
        'id': 'GLC_total_input',
        'comp': [('GLC_input', ), ('GLC_unlabelled_input', )],
        ModelKeyword.flux_range: ModelKeyword.add_range_type
    }
]

