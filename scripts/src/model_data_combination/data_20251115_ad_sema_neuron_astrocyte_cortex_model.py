from .data_20251115_ad_sema import (
    DataSource, common_data_loader, model_loader, ModelList, common_model_constructor,
    Keywords, MFARunningMode, it, ModelKeyword,
    UserDefinedMFASettings as OldUserDefinedMFASettings, MFASettings,
    result_process_information_dict_analysis, flux_process_map_func,
    result_process_figure_config_dict_generator, display_flux_name_dict as old_display_flux_name_dict,
    metabolic_network_parameter_generator, dataset_extractor,
)

data_name = DataSource.data_20251115_ad_sema
model_name = ModelList.liver_neuron_astrocyte_cortex_model
_, keyword = common_data_loader(
    data_name, keyword_only=True, test_mode=False, natural_anti_correction=False)
(
    user_defined_model, labeled_buffer_input_metabolite_obj_dict,
    specific_reaction_range_dict) = model_loader(model_name)
mfa_model_obj = common_model_constructor(user_defined_model)

original_tissue_name_tuple = (ModelKeyword.cortex,)
new_cortex_tissue_name_tuple = (ModelKeyword.astrocyte, ModelKeyword.neuron)
def return_data_content(running_mode, dataset_id=None, **kwargs):
    data_content, _ = common_data_loader(data_name, test_mode=False, natural_anti_correction=False)
    if running_mode == MFARunningMode.raw_experimental_data_plotting:
        return data_content
    else:
        (mfa_data_dict, data_content_1, data_content_2, mfa_information_dict, data_content_4) = data_content
        complete_dataset_dict = {
            'part1': [
                'glucose__ad__ctrl__1',
                # 'glucose__ad__ctrl__2',
                # 'glucose__ad__ctrl__3',
                'glucose__ad__semaglutide__1',
                # 'glucose__ad__semaglutide__2',
                # 'glucose__ad__semaglutide__3',
            ],
            'part2': [
                # 'glucose__ad__semaglutide__1',
                'glucose__ad__ctrl__2',
                'glucose__ad__ctrl__3',
                'glucose__ad__ctrl__4',
                'glucose__ad__ctrl__5',
                'glucose__ad__ctrl__6',
            ],
            'part3': [
                'glucose__ad__semaglutide__2',
                'glucose__ad__semaglutide__3',
                'glucose__ad__semaglutide__4',
                'glucose__ad__semaglutide__5',
                'glucose__ad__semaglutide__6',
                # 'glucose__ad__ctrl__6',
            ]
        }
        new_mfa_data_dict, new_mfa_information_dict = dataset_extractor(
            mfa_data_dict, mfa_information_dict, dataset_id, complete_dataset_dict)
        for mfa_data_obj in new_mfa_data_dict.values():
            mfa_data_obj.add_input_metabolite_obj_dict(labeled_buffer_input_metabolite_obj_dict)
            mfa_data_obj.replace_tissue_in_mid_data(original_tissue_name_tuple, new_cortex_tissue_name_tuple, copy=False)
        new_data_content = (new_mfa_data_dict, data_content_1, data_content_2, new_mfa_information_dict, data_content_4)
        return new_data_content


class UserDefinedMFASettings(OldUserDefinedMFASettings):
    tf2_solver = True
    optimization_num_for_each_suffix_dict = {
        None: 10000,
    }
    user_defined_constant_flux_dict = {
        **MFASettings.user_defined_constant_flux_dict,
        'GLC_infuse': (ModelKeyword.serum, 200)
    }
    user_defined_specific_flux_range_dict = {
        **MFASettings.user_defined_specific_flux_range_dict,
        'Salvage_c': (None, (1, 10)),
        'PROS_m': (None, (1, 5)),

        'BHB_input': (ModelKeyword.astrocyte, (1, 20)),
        'SCOT_m': (ModelKeyword.astrocyte, (1, 10)),
        'ACSS2_m': (ModelKeyword.astrocyte, (1, 50)),
        'ACAC_production_m': (ModelKeyword.astrocyte, (5, 200)),
        'PC_m': (ModelKeyword.astrocyte, (5, 200)),
        'GLNS_c': (ModelKeyword.astrocyte, (10, 300)),
        'FA_usage': (ModelKeyword.astrocyte, (1, 600)),
        'GLND_m': (ModelKeyword.astrocyte, (1, 100)),
        'ASP_exchange': (ModelKeyword.astrocyte, (1, 10)),
        'ASP_exchange__R': (ModelKeyword.astrocyte, (1, 10)),

        'ALA_c_isf_exchange__R': {ModelKeyword.neuron: (1, 50), ModelKeyword.astrocyte: (1, 50)},
        'ASP_c_isf_exchange': {ModelKeyword.neuron: (1, 50), ModelKeyword.astrocyte: (1, 50)},
        'ASP_c_isf_exchange__R': {ModelKeyword.neuron: (1, 50), ModelKeyword.astrocyte: (1, 50)},
        'GLYCEROL_exchange': (ModelKeyword.astrocyte, (1, 100)),
        'GLYCEROL_exchange__R': (ModelKeyword.astrocyte, (1, 100)),
        'HEX_c': {ModelKeyword.neuron: (70, 600), ModelKeyword.astrocyte: (10, 600)},
        'CS_m': {ModelKeyword.neuron: (100, 1000), ModelKeyword.astrocyte: (50, 900)},
        'PDH_m': {ModelKeyword.neuron: (100, 800), ModelKeyword.astrocyte: (10, 600)},
        'BHBDH_m': {ModelKeyword.neuron: (1, 50), ModelKeyword.astrocyte: (1, 50)},
        'BHBDH_m__R': {ModelKeyword.neuron: (1, 50), ModelKeyword.astrocyte: (1, 50)},
        'BHB_trans': {ModelKeyword.neuron: (1, 100), ModelKeyword.astrocyte: (1, 100)},
        'BHB_trans__R': {ModelKeyword.neuron: (1, 100), ModelKeyword.astrocyte: (1, 100)},
        'AKGMAL_m': {ModelKeyword.neuron: (1, 1000), ModelKeyword.astrocyte: (1, 1000)},
        'AKGMAL_m__R': {ModelKeyword.neuron: (1, 1000), ModelKeyword.astrocyte: (1, 1000)},
        'ME2_c': {ModelKeyword.neuron: (1, 30), ModelKeyword.astrocyte: (1, 30)},

        'ACOT_c': (ModelKeyword.neuron, (1, 100)),
        'PEPCK_c': (ModelKeyword.neuron, (1, 100)),
        'LAC_c_isf_exchange__R': (ModelKeyword.neuron, (1, 600)),
        'ASPGLU_m': (ModelKeyword.neuron, (1, 1000)),
        'ASPGLU_m__R': (ModelKeyword.neuron, (1, 1000)),
    }
    user_defined_specific_flux_range_dict = MFASettings.combine_flux_dict_list(
        user_defined_specific_flux_range_dict, specific_reaction_range_dict.values())
    average_optimized_results = None

    def tf2_running_parameters(self):
        batch_size_32g = 120
        batch_size_64g = 150
        batch_size_96g = 240
        nnls_parallel_num_32g_64g = 30
        nnls_parallel_num_96g = 40
        emu_gh_parallel_num_default = 3
        test_batch_size = 39
        parameter_tuple_for_each_suffix_dict = {
            'test': (test_batch_size, 2000, test_batch_size, test_batch_size, emu_gh_parallel_num_default),
            None: (batch_size_32g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g, emu_gh_parallel_num_default),
            '0': (batch_size_32g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g, emu_gh_parallel_num_default),
            '1': (batch_size_64g, 4000, nnls_parallel_num_32g_64g, nnls_parallel_num_32g_64g, emu_gh_parallel_num_default),
            '2': (batch_size_96g, 4000, nnls_parallel_num_96g, nnls_parallel_num_96g, emu_gh_parallel_num_default),
        }
        return batch_size_32g, parameter_tuple_for_each_suffix_dict

user_defined_mfa_settings = UserDefinedMFASettings()

display_flux_name_dict = {
    **old_display_flux_name_dict,
    ('LAC_c_isf_exchange', 'LAC_c_isf_exchange__R'): 'LAC isf output net',
    ('PYR_c_isf_exchange', 'PYR_c_isf_exchange__R'): 'PYR isf output net',
    ('ALA_c_isf_exchange', 'ALA_c_isf_exchange__R'): 'ALA isf output net',
    'GLU_c_isf_exchange': 'GLU isf output net',
    'GLN_c_isf_exchange__R': 'GLN isf input net',
}



def tissue_mid_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    mid_name_list_in_serum = [
        ['glucose', 'pyruvate', 'lactate', 'alanine', 'serine', ],
        ['glutamate', 'glutamine', 'aspartate', '3-hydroxybutyrate'],
    ]
    mid_name_list_in_one_tissue = [
        [
            'glucose', 'pyruvate', 'lactate', 'alanine', 'acetyl-coa', '3-hydroxybutyrate'],
        [
            'citrate', 'glutamate', 'glutamine', 'succinate', 'fumarate',
            'malate', 'aspartate', ],
    ]

    base_tissue_mid_name_list_dict = {
        ModelKeyword.serum: mid_name_list_in_serum,
        ModelKeyword.liver: mid_name_list_in_one_tissue,
        new_cortex_tissue_name_tuple: mid_name_list_in_one_tissue,
    }
    tissue_mid_name_list_dict = base_tissue_mid_name_list_dict
    return tissue_mid_name_list_dict, keyword.display_mid_name_dict, keyword.display_tissue_name_dict


def tissue_flux_name_list_dict_constructor(current_model_name=model_name, **kwargs):
    flux_name_list_in_serum = [
        ['GLYCOLYSIS_c', ('LDH_e', 'LDH_e__R'), ],
    ]

    flux_name_list_in_liver = [
        [('GLC_exchange', 'GLC_exchange__R'), ('FBA_c', 'FBA_c__R'), 'PYK_c', ('LDH_c', 'LDH_c__R'),
         ('LAC_exchange', 'LAC_exchange__R'), 'BHB_output'],
        ['CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), ('ALA_exchange', 'ALA_exchange__R')],
    ]

    flux_name_list_in_astrocyte = [
        [('GLC_exchange', 'GLC_exchange__R'), ('FBA_c', 'FBA_c__R'), 'PYK_c', ('LDH_c', 'LDH_c__R'),
         ('LAC_exchange', 'LAC_exchange__R'), 'PHGDH_PSAT_PSP_c', ('SER_exchange', 'SER_exchange__R')],
        [('GLYCEROL_exchange', 'GLYCEROL_exchange__R'), 'BHB_input', 'FA_usage',
         'CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), ('ALA_exchange', 'ALA_exchange__R')],
    ]

    flux_name_list_in_neuron = [
        [('GLC_exchange', 'GLC_exchange__R'), ('FBA_c', 'FBA_c__R'), 'PYK_c', ('LDH_c', 'LDH_c__R'),
         ('LAC_c_isf_exchange', 'LAC_c_isf_exchange__R'), 'PDH_m', 'SER_c_isf_exchange__R', ('ACSS2_m', 'ACOT_c'), 'FA_usage'],
        ['GLU_c_isf_exchange', 'GLN_c_isf_exchange__R', ('PYR_c_isf_exchange', 'PYR_c_isf_exchange__R'),
         'CS_m', 'AKGD_m', ('MDH_m', 'MDH_m__R'), ('ALA_c_isf_exchange', 'ALA_c_isf_exchange__R'), ('ACE_trans', 'ACE_trans__R')],
    ]

    all_tissue_flux_name_list_dict = {
        ModelKeyword.serum: flux_name_list_in_serum,
        ModelKeyword.liver: flux_name_list_in_liver,
        ModelKeyword.astrocyte: flux_name_list_in_astrocyte,
        ModelKeyword.neuron: flux_name_list_in_neuron,
    }
    specific_tissue_reversible_flux_pair_dict = {
    }
    return (
        all_tissue_flux_name_list_dict, display_flux_name_dict, keyword.display_tissue_name_dict,
        specific_tissue_reversible_flux_pair_dict, flux_process_map_func)
