ketone_body_utilization_reaction_list = [
    {
        'id': 'BHB_input',
        'sub': [('BHB_e', 'abcd')],
        'pro': [('BHB_c', 'abcd')],
    },
    {
        'id': 'BHB_trans',
        'sub': [('BHB_c', 'abcd')],
        'pro': [('BHB_m', 'abcd')],
        'reverse': True
    },
    {
        'id': 'BHBDH_m',
        'sub': [('BHB_m', 'abcd')],
        'pro': [('ACAC_m', 'abcd')],
        'reverse': True
    },
    {
        'id': 'SCOT_m',
        'sub': [('ACAC_m', 'abcd')],
        'pro': [('ACACCOA_m', 'abcd')],
    },
    {
        'id': 'ACACCOA_lyase_m',
        'sub': [('ACACCOA_m', 'abcd')],
        'pro': [('ACCOA_m', 'ab'), ('ACCOA_m', 'cd')],
    },
]

extra_reaction_list = [
    *ketone_body_utilization_reaction_list,
]

emu_excluded_metabolite_set = set()

balance_excluded_metabolite_set = emu_excluded_metabolite_set

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set,
    ModelKeyword.serum: {'BHB_e', }
}

