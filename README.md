Code used to generate results for the paper: *Integrated DNA methylation and gene expression profiling across multiple brain regions implicate novel genes in Alzheimer’s disease*

## Citation
Coming soon

## Script summary

#### Importing idats, preprocessing, and exploratory data analysis
* *00_load_idats.R*:
* 01a_preprocess_methylation.R
* 01b_check_genotype_correlation_meth450k_vs_SNPChip.R
* 01c_check_negative_control_PCs.R
* 01d_pca.R
* 01e_subset_to_samples_for_analysis.R
* 01f_pull_genotypes.R
* 01g_create_demographic_table.R
* 01h_horvath_epigenetic_clock.R

#### Probe-level analyses, results summary, and visualizations
* 02a_i_caseControl_regionalSpecific_probe_differences.R
* 02a_iii_caseControl_mainResultsModel_DMP_NeuN_sensitivity.R
* 02a_iv_caseControl_mainResultsModel_DMP_APOE_sensitivity.R
* 02a_v_normal_aging_controlsOnly.R
* 02a_vi_caseControl_noCRB_crossRegion_DMP.R
* 02b_merge_allRegion_stats.R
* 02c_compare_DMP_results_to_postmortem_AD_studies.R
* 02d_caseControl_ageAcceleration_and_Composition.R
* 02e_caseControl_boxplots.R
* 02f_caseControl_gene_set_enrichment.R
* 02g_caseControl_heatmaps_DMP.R
* 02h_sensitivity_posthoc_distribution_plots.R

#### Region-level analyses, results summary, and visualizations
* 03a_caseControl_DMR_analysis.R
* 03b_caseControl_DMR_plots.R
* 03c_caseControl_DMR_analysis_NeuN_sensitivity.R

#### Integrating gene-expression data and asessing replicability
* 04a_comparing_global_alz_DMP_stats_across_datasets.R
* 04b_check_case_control_stats_for_DMP_genes.R
* 04c_correlate_DNAm_with_gene_expression.R
* 04d_replicability_of_top_DMPs.R

#### Checking involvement of aging
* 05a_aging_control_AD_Int.R
* 05c_aging_results_analysis.R

#### Functional and system-level analyses
* 06a_genetic_risk_loci_DMP_results.R
* get_cpg_in_risk_loci.R
* scrape_nature_genetics_lambert_et_al.R
* 06b_string_ppi_networks.R
* 06c_check_coexpression_networks.R

#### Reprocessing data from Lunnon et al. (2014).
* Lunnon_2014_getGEO_01.R
* Lunnon_2014_DMP_02.R
