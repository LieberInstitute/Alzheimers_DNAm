################################# Run all analysis with LD block loci
ad_genetic_risk_loci = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_LD_block_loci_cpgs.csv')
load('rdas/caseControl_DMC_allRegion.rda')
cpg_in_loci = cbind(allStats[as.character(ad_genetic_risk_loci$Name),],ad_genetic_risk_loci[,-1] )
table(is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val))
cpg_in_loci = cpg_in_loci[!is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val),]

## Create summary table
library(dplyr)
library(tidyr)
cpg_in_loci = group_by(cpg_in_loci,indexSnp_nearestGene )
genetic_risk_DMP_summary = summarize(cpg_in_loci, 
									 GeneticMarker = unique(indexSnp_rsNum),
									 Chromosome = unique(chr_hg38), 
									 Covered_Region = paste0(min(pos_hg38),":",max(pos_hg38)) , 
									 Region_Length = max(pos_hg38) - min(pos_hg38), 
									 Number_Of_CpGs = n(), 
									 n_mainEffect_nomP05= sum(Primary_subset_mainEffect_P.Value<0.05),
									 n_mainEffect_FDR05= sum(Primary_subset_mainEffect_adj.P.Val<0.05),
									 n_interaction_nomP05=sum(Primary_subset_interactionEffect_P.Value<0.05),
									 n_interactionEffect_FDR05= sum(Primary_subset_interactionEffect_adj.P.Val<0.05)) %>%
									 as.data.frame()
write.csv(genetic_risk_DMP_summary, file='csvs/SupplementalTable_genetic_risk_dmp_overlap_summary_ld_block_loci.csv',row.names=F )

### Check enrichment
allStats$within_genetic_risk_loci = allStats$Name%in% ad_genetic_risk_loci$Name
######## Main effect
### All
fisher.test(table(riskLoci=allStats$within_genetic_risk_loci,DMP=allStats$Primary_subset_mainEffect_adj.P.Val<0.05 ) ) ## FDR cutoff
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_mainEffect_P.Value<0.01 ) )
### hypermethylated
hyper = allStats[allStats$Primary_subset_mainEffect_logFC>0,]
fisher.test(table(riskLoci=hyper$within_genetic_risk_loci,DMP=hyper$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(hyper$within_genetic_risk_loci,hyper$Primary_subset_mainEffect_P.Value<0.05 ) )

### hypomethylated
hypo = allStats[allStats$Primary_subset_mainEffect_logFC<0,]
fisher.test(table(riskLoci=hypo$within_genetic_risk_loci,DMP=hypo$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(hypo$within_genetic_risk_loci,hypo$Primary_subset_mainEffect_P.Value<0.05 ) )

######## interaction effect
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_P.Value<0.01 ) )

##
fisher.test(table(riskLoci=allStats$within_genetic_risk_loci,DMP=allStats$Sensitivity_subset_mainEffect_adj.P.Val<0.05 ) ) ## FDR cutoff
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Sensitivity_subset_mainEffect_P.Value<0.01 ) )

fisher.test(table(allStats$within_genetic_risk_loci,allStats$Sensitivity_subset_interactionEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Sensitivity_subset_interactionEffect_P.Value<0.01 ) )
