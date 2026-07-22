from .common.packages import ValueEnum


class ExperimentName(ValueEnum):
    data_20251115_ad_sema = 'data_20251115_ad_sema'
    data_20251115_ad_sema_serum_labeled_buffer = 'data_20251115_ad_sema_serum_labeled_buffer'
    data_20251115_ad_sema_single_tissue_model = 'data_20251115_ad_sema_single_tissue_model'
    data_20251115_ad_sema_serum_labeled_buffer_tf2 = 'data_20251115_ad_sema_serum_labeled_buffer_tf2'
    data_20251115_ad_sema_neuron_astrocyte_cortex_model = 'data_20251115_ad_sema_neuron_astrocyte_cortex_model'

# default_data_model_name = ExperimentName.data_20251115_ad_sema
# default_data_model_name = ExperimentName.data_20251115_ad_sema_serum_labeled_buffer
default_data_model_name = ExperimentName.data_20251115_ad_sema_single_tissue_model
# default_data_model_name = ExperimentName.data_20251115_ad_sema_serum_labeled_buffer_tf2
# default_data_model_name = ExperimentName.data_20251115_ad_sema_neuron_astrocyte_cortex_model


class MetabolomicsExperimentName(ValueEnum):
    data_20260412_metabolomics = 'data_20260412_metabolomics'


class MetabolomicsRunningMode(ValueEnum):
    summary = 'summary'
    export = 'export'
    enrichment = 'enrichment'
    plsda = 'plsda'
    oplsda = 'oplsda'
    clustered_heatmap = 'clustered_heatmap'
    pathway_bubble = 'pathway_bubble'
    splot = 'splot'
    phenotype_corr = 'phenotype_corr'
    volcano = 'volcano'


class MLRunningMode(ValueEnum):
    met_gat = 'met_gat'


class MFARunningMode(ValueEnum):
    flux_analysis = 'flux_analysis'
    raw_experimental_data_plotting = 'raw_experimental_data_plotting'
    result_process = 'result_process'
    solver_output = 'solver_output'
    flux_analysis_test = 'flux_analysis_test'


data_model_comment = {
    ExperimentName.data_20251115_ad_sema:
        'Data from mice',
    ExperimentName.data_20251115_ad_sema_serum_labeled_buffer:
        'Data from mice',
    ExperimentName.data_20251115_ad_sema_single_tissue_model:
        'Data from mice',
    ExperimentName.data_20251115_ad_sema_serum_labeled_buffer_tf2:
        'Data from mice',
    ExperimentName.data_20251115_ad_sema_neuron_astrocyte_cortex_model:
        'Data from mice',
}

metabolomics_data_model_comment = {
    MetabolomicsExperimentName.data_20260412_metabolomics:
        'AD + Semaglutide multi-tissue metabolomics',
}
