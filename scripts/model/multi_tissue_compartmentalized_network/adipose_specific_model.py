extra_reaction_list = [
    {
        # The carbon atom mapping is just for simplification of emu analyzer
        'id': 'lipid_mobilization_c',
        'sub': [('LIPID', 'abcde')],
        'pro': [('GLYCEROL_c', 'abc'), ('FA_c', 'de', 12)],
    },
]

emu_excluded_metabolite_set = {'LIPID'}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set
}

