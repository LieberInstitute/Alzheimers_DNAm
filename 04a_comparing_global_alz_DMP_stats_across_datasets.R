## Comparing global alzheimer statistics 
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
## Load our statistics
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')

load('rdas/tidyStats_caseControl_DMC_singleRegion.rda')
load('rdas/caseControl_DMC_allRegion.rda',verbose=T)

## Load watson et al 2016: STG
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE76105/GSE76105_Watson_mergedStats.rda')
colnames(Watson_mergedStats) = paste0(colnames(Watson_mergedStats), "_", "Watson2016")

## Load lunnon et al 2014: Multiregion
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_Lunnon_mergedStats.rda')
colnames(Lunnon_mergedStats) = paste0(colnames(Lunnon_mergedStats), "_", "Lunnon2014")

## Load Smith et al 2018
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE80970/GSE80970_Smith_mergedStats.rda',verbose=T)
colnames(Smith_mergedStats) = paste0(colnames(Smith_mergedStats), "_", "Smith2018")

##
ii = rownames(allStats)
multisetStats = cbind(allStats,
					  Watson_mergedStats[match(ii, rownames(Watson_mergedStats)),],
					  Lunnon_mergedStats[match(ii, rownames(Lunnon_mergedStats)),],
					  Smith_mergedStats[match(ii, rownames(Smith_mergedStats)),])


### Checking enrichments

## Lunnon et al 2014
lun = table(multisetStats$Primary_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$Primary_subset_mainEffect_logFC) )
fisher.test(lun)

## Smith et al 2018
smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,
multisetStats$ALL_main_Braak_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Braak_logFC_Smith2018) == sign(multisetStats$ALL_subset_mainEffect_logFC) )
fisher.test(smith)

smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Dx_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) )

smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Braak_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Braak_logFC_Smith2018) == sign(multisetStats$ALL_subset_mainEffect_logFC) )

## Watson et al 2016
watson = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$STG_Dx_Unadjusted_P.Value_Watson2016<0.05 & sign(multisetStats$STG_Dx_Unadjusted_logFC_Watson2016) == sign(multisetStats$ALL_subset_mainEffect_logFC) )

###
anySet = watson = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05, 
(multisetStats$STG_Dx_Unadjusted_P.Value_Watson2016<0.05 & sign(multisetStats$STG_Dx_Unadjusted_logFC_Watson2016) == sign(multisetStats$ALL_subset_mainEffect_logFC) ) | 
(multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC)) | 
(multisetStats$ALL_main_Dx_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) ) ) 

topHits = multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & 
(multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC))  &multisetStats$ALL_main_Dx_P.Value_Smith2018<0.05& sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) 
multisetStats[topHits,'within10kb_geneSymbol_gencode_hg38']
####### Saving out replication results

write.csv(multisetStats[multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC),],file='csvs/SupplementalTable_dmp_replicated_nomSig_dirConsistent_in_lunnon2014.csv',row.names=F )

###
j = table(multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05,multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025  )
fisher.test(j)

write.csv(multisetStats[multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05&multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025,], file='csvs/SupplementalTable_dmp_replicated_nomSig_interaction_in_lunnon2014.csv',row.names=F )
