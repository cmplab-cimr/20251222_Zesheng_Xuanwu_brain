from scripts.src.core.common.functions import exchange_reaction_constructor
from scripts.src.core.common.config import ModelKeyword

extra_ketone_body_production_reaction_list_for_astrocyte = [
    {
        'id': 'ACAC_production_m',
        'sub': [('ACCOA_m', 'ab'), ('ACCOA_m', 'cd')],
        'pro': [('ACAC_m', 'abcd')],
    },
]

isf_exchange_metabolite_carbon_num_list = [
    ('GLU_c', 5),
    ('GLN_c', 5),
    ('SER_c', 3),
    ('LAC_c', 3),
    ('ALA_c', 3),
    ('PYR_c', 3),
    ('ASP_c', 4),
    ('BHB_c', 4),
    # ('GLC_e', 6),
    # ('GLY_e', 2),
    # ('ACE_e', 2),
    # ('BHB_e', 4),
    # ('GLYCEROL_e', 3),
    # ('FA_e', 2),
]

lowered_flux_range = (1, 5)
neuron_delete_reaction_dict = {
    # Cellular
    'PC_m': 'F',                   # 丙酮酸羧化酶:星形胶质细胞特异,神经元几乎不表达 -> 关闭
    'GLNS_c': 'F',                 # 谷氨酰胺合成酶:星形胶质细胞特异 -> 关闭
    # 'PHGDH_PSAT_PSP_c': 'F',       # de novo 丝氨酸合成(PHGDH):星形胶质细胞特异;神经元改为从 ISF 摄取丝氨酸 -> 关闭

    # Serum
    # GLC_exchange 开放(GLUT3 直接摄取血糖)。
    'GLU_exchange': 'B',   # 谷氨酸不从血浆喂入(BBB 几乎不通)
    'GLN_exchange': 'B',   # 谷氨酰胺从 ISF 摄取(星胶供给),不直接取血浆
    'SER_exchange': 'B',   # 丝氨酸从 ISF 摄取(星胶 PHGDH 供给)
    'LAC_exchange': 'B',   # 乳酸走 ISF(ANLS),不直接碰血浆
    'PYR_exchange': 'B',   # 丙酮酸走 ISF
    'ALA_exchange': 'B',   # 丙氨酸走 ISF
    'ASP_exchange': 'B',   # 天冬氨酸非血源底物
    'ACE_exchange': 'B',   # 醋酸是胶质特异底物,神经元不摄取
    'BHB_input': 'F',      # 神经元不直接从血液中摄取酮体

    # ISF
    'GLU_c_isf_exchange': 'R',      # 神经元净【释放】谷氨酸
    'GLN_c_isf_exchange': 'F',      # 净【摄取】谷氨酰胺/丝氨酸
    'SER_c_isf_exchange': 'F',
}
neuron_specific_reaction_range_dict = {
    key: {ModelKeyword.neuron: lowered_flux_range}
    for key in [
        # 'PC_m',                   # 丙酮酸羧化酶:星形胶质细胞特异,神经元几乎不表达 -> 关闭
        # 'GLNS_c',                 # 谷氨酰胺合成酶:星形胶质细胞特异 -> 关闭
        'PHGDH_PSAT_PSP_c',       # de novo 丝氨酸合成(PHGDH):星形胶质细胞特异;神经元改为从 ISF 摄取丝氨酸 -> 关闭
        'GLYCEROL_exchange',      # 甘油/脂肪酸:神经元主要用于膜脂而非供能,这里按"碳源=葡萄糖"关闭
        'GLYCEROL_exchange__R',   # 甘油/脂肪酸:同上,可逆反应的反向
        'FA_exchange',            # (若要建模神经元脂质摄取可重新打开这两条)
        'FA_exchange__R',
    ]
}
astrocyte_delete_reaction_dict = {
    # Cellular
    'ASPGLU_m': 'B',  # 天冬氨酸转氨酶(线粒体)

    # Serum
    'GLU_exchange': 'B',  # 谷氨酸不进/出血浆;只在脑内循环(神经元->ISF->星胶)

    # ISF
    'GLU_c_isf_exchange': 'F',     # 星胶净【摄取】谷氨酸(EAAT/GLUt6)-> 关释放臂
    'GLN_c_isf_exchange': 'R',     # 净【释放】谷氨酰胺/丝氨酸 -> 关再摄取臂
    'SER_c_isf_exchange': 'R',
}
astrocyte_specific_reaction_range_dict = {
    key: {ModelKeyword.astrocyte: lowered_flux_range}
    for key in [
        # 'ASPGLU_m',     # 天冬氨酸转氨酶(线粒体)
        # 'ASPGLU_m__R',  # 同上,可逆反应的反向
    ]
}

isf_exchange_reaction_list, isf_metabolite_list = exchange_reaction_constructor(
    isf_exchange_metabolite_carbon_num_list, original_compartment='c', target_compartment=ModelKeyword.isf,
    reaction_specific_prefix=None, reaction_specific_suffix='exchange', reversible=True)

extra_reaction_list = isf_exchange_reaction_list

emu_excluded_metabolite_set = {'BIOMASS',}

balance_excluded_metabolite_set = emu_excluded_metabolite_set

symmetrical_metabolite_set = set()

model_compartment_set = {'e', ModelKeyword.isf}

composite_reaction_list = []

other_pool_metabolite_dict = {
    ModelKeyword.general: emu_excluded_metabolite_set,
    ModelKeyword.isf: set(isf_metabolite_list),
}

