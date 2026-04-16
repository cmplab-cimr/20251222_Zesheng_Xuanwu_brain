from .liver_specific_compartmentalized_model import (
extra_reaction_list as liver_extra_reaction_list,
emu_excluded_metabolite_set,
balance_excluded_metabolite_set,
other_pool_metabolite_dict
)

extra_reaction_list = [
    reaction_dict for reaction_dict in liver_extra_reaction_list
    if reaction_dict['id'] != 'GLYCOGENOLYSIS_c'
]
