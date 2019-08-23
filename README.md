Code used to generate results for the paper: *Integrated DNA methylation and gene expression profiling across multiple brain regions implicate novel genes in Alzheimer's disease*

[![DOI](https://zenodo.org/badge/133691598.svg)](https://zenodo.org/badge/latestdoi/133691598)

## Citation
Semick SA, Bharadwaj RA, Collado-Torres L, Tao R, Shin JH, Deep-Soboslay A, Weiss JR, Weinberger DR, Hyde TM, Kleinman JE, Jaffe AE, Mattay VS. Integrated DNA methylation and gene expression profiling across multiple brain regions implicate novel genes in Alzheimer's disease. 2019 Feb 2. PMID: 30712078 

URL: https://link.springer.com/article/10.1007%2Fs00401-019-01966-5


## Script summary

#### Importing idats, preprocessing, and exploratory data analysis
* *00_load_idats.R*: Import idat files into minfi (an RGset).
* *01a_preprocess_methylation.R*: Preprocess methylation array data (drop low-quality samples and probes, check QC).
* *01b_check_genotype_correlation_meth450k_vs_SNPChip.R*: Check correlation between SNP-chip and methylation array genotypes (~65 probes).
* *01c_check_negative_control_PCs.R*: Inspect principal component analysis of negative control probes.
* *01d_pca.R*: Principal component analysis of probes used in analysis.
* *01e_subset_to_samples_for_analysis.R*: Removal of samples not used in analysis.
* *01f_pull_genotypes.R*: [Not used in paper].
* *01g_create_demographic_table.R*: Create a demographic table (Supplemental Table 1).
* *01h_horvath_epigenetic_clock.R*: Compute the estimated DNAm age via Horvath's clock.

#### Probe-level analyses, results summary, and visualizations
* *02a_i_caseControl_regionalSpecific_probe_differences.R*: Case-control differential methylated probe (DMP) analysis stratified by brain region.
* *02a_iii_caseControl_mainResultsModel_DMP_NeuN_sensitivity.R*: Case-control DMP analysis, adjusting for estimates of neuronal proportion.
* *02a_iv_caseControl_mainResultsModel_DMP_APOE_sensitivity.R*: Case-control DMP analysis, adjusting APOE4 dosage.
* *02a_v_normal_aging_controlsOnly.R*: DMP effect of aging on DNAm in normal controls.
* *02a_vi_caseControl_noCRB_crossRegion_DMP.R*: [Not used in paper].
* *02b_merge_allRegion_stats.R*: Merging together statistics from multiple models.
* *02c_compare_DMP_results_to_postmortem_AD_studies.R*: Assessing replication of previously reported DMPs in our dataset.
* *02d_caseControl_ageAcceleration_and_Composition.R*: Case-control differences for DNAm age acceleration and for neuronal cell type proportions
* *02e_caseControl_boxplots.R*: Boxplots of top DMPs via various models.
* *02f_caseControl_gene_set_enrichment.R*: Enrichment of DNAm probes in genes.
* *02g_caseControl_heatmaps_DMP.R*: Heatmaps of top DMPs.
* *02h_sensitivity_posthoc_distribution_plots.R*: Posthoc distribution sensitivity plots.

#### Region-level analyses, results summary, and visualizations
* *03a_caseControl_DMR_analysis.R*: Run DNAm region-level case-control analysis.
* *03b_caseControl_DMR_plots.R*: Plot significant differentially methylated regions.
* *03c_caseControl_DMR_analysis_NeuN_sensitivity.R*: Run DMR case-control analysis with NeuN estimates adjustment.

#### Integrating gene-expression data and asessing replicability
* *04a_comparing_global_alz_DMP_stats_across_datasets.R*: Check correlation between statistics across datasets.
* *04b_check_case_control_stats_for_DMP_genes.R*: Pull differential gene expression results for case-control analysis.
* *04c_correlate_DNAm_with_gene_expression.R*: Correlate DNAm with gene expression of nearby genes (within 10kb either side).
* *04d_replicability_of_top_DMPs.R*: Check if top DMPs are present in Lunnon et al. 2014.

#### Checking involvement of aging
* *05a_aging_control_AD_Int.R*: Not included in paper.
* *05c_aging_results_analysis.R*: Comparing "normal" (unaffected-control) aging  DNAm changes to AD-associated DNAm changes.

#### Functional and system-level analyses
* *06a_genetic_risk_loci_DMP_results.R*: Assess enrichment of DMPs in GWAS loci.
* *get_cpg_in_risk_loci.R*: Determine which CpG probes lie within GWAS risk loci.
* *scrape_nature_genetics_lambert_et_al.R*: Pull index SNP from 2014 AD-GWAS (Lambert et al.) to determine AD risk loci. 
* *06b_string_ppi_networks.R*: Not included in paper.
* *06c_check_coexpression_networks.R*: Not included in paper.

#### Reprocessing data from Lunnon et al. (2014).
* *Lunnon_2014_getGEO_01.R*: Download Lunnon et al data.
* *Lunnon_2014_DMP_02.R*: Model case-control differences.
