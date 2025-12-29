extra_reaction_list = [
    {
        'id': 'PGMT_c',
        'sub': [('GLC6P_c', 'abcdef')],
        'pro': [('GLC1P_c', 'abcdef')],
        'reverse': True
    },
    {
        'id': 'GLYCOGENOLYSIS_c',
        'sub': [('GLYGN_e', 'abcdef')],
        'pro': [('GLC1P_c', 'abcdef')],
    },
]

emu_excluded_metabolite_set = {'GLYGN_e'}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set
}

