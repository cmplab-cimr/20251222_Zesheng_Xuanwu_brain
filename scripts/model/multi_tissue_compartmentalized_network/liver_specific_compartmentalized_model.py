extra_glycogen_reaction_list = [
    {
        'id': 'G6Pase_c',
        'sub': [('GLC6P_c', 'abcdef')],
        'pro': [('GLC_c', 'abcdef')],
    },
    {
        'id': 'FBPase_c',
        'sub': [('FRU16BP_c', 'abcdef')],
        'pro': [('FRU6P_c', 'abcdef')],
    },
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


extra_fatty_acid_reaction_list = [
]

extra_ketone_body_reaction_list = [
    {
        'id': 'ACAC_production_m',
        'sub': [('ACCOA_m', 'ab'), ('ACCOA_m', 'cd')],
        'pro': [('ACAC_m', 'abcd')],
    },
    {
        'id': 'BHBDH_m',
        'sub': [('BHB_m', 'abcd')],
        'pro': [('ACAC_m', 'abcd')],
        'reverse': True
    },
    {
        'id': 'BHB_trans',
        'sub': [('BHB_c', 'abcd')],
        'pro': [('BHB_m', 'abcd')],
        'reverse': True
    },
    {
        'id': 'BHB_output',
        'sub': [('BHB_c', 'abcd')],
        'pro': [('BHB_e', 'abcd')],
    },
]
extra_reaction_list = [
    *extra_glycogen_reaction_list,
    *extra_fatty_acid_reaction_list,
    *extra_ketone_body_reaction_list,
]

emu_excluded_metabolite_set = {'GLYGN_e', 'FA_PRP_c', }

balance_excluded_metabolite_set = emu_excluded_metabolite_set

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set,
    ModelKeyword.serum: {'BHB_e', 'ACE_e', }
}

