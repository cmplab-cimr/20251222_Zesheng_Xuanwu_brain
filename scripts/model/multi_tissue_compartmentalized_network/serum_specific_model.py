reaction_list = [
    {
        'id': 'GLC_infuse',
        'sub': [('GLC_inf', 'abcdef')],
        'pro': [('GLC_e', 'abcdef')],
    },
    {
        'id': 'GLYCOLYSIS_c',
        'sub': [('GLC_e', 'abcdef')],
        'pro': [('PYR_e', 'cba'), ('PYR_e', 'def')],
    },
    {
        'id': 'LDH_e',
        'sub': [('PYR_e', 'abc')],
        'pro': [('LAC_e', 'abc')],
        'reverse': True
    },
]

emu_excluded_metabolite_set = {'GLC_inf'}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

added_input_metabolite_set = {'GLC_inf'}

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set
}
