
from .basic_tissue_compartmentalized_model import (
    other_pool_metabolite_dict, ModelKeyword, biomass_reaction_dict_constructor)
from scripts.src.core.common.functions import exchange_reaction_constructor


serum_extra_metabolite_set = {key for key in other_pool_metabolite_dict[ModelKeyword.serum] if key not in {'GLY_e'}}

extra_reaction_list, buffer_metabolite_list = exchange_reaction_constructor(
    [
        ('GLC_e', 6),
        ('GLU_e', 5),
        ('GLN_e', 5),
        ('ASP_e', 4),
        ('SER_e', 3),
        # ('GLY_e', 2),
        ('ALA_e', 3),
        ('LAC_e', 3),
        ('PYR_e', 3),
        ('ACE_e', 2),
        ('BHB_e', 4),
        ('GLYCEROL_e', 3),
        ('FA_e', 2),
    ], original_compartment='e', target_compartment='buffer',
    reaction_specific_suffix=None, reversible=True)

emu_excluded_metabolite_set = {'BIOMASS', *buffer_metabolite_list}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

symmetrical_metabolite_set = set()

model_compartment_set = {'e'}

composite_reaction_list = []

from scripts.src.core.common.config import ModelKeyword
other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set,
    # ModelKeyword.serum: serum_extra_metabolite_set,
}

