library(RColorBrewer)
library(pheatmap)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tidyr)
library(dplyr)
library(data.table)

col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
## 
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')

primaryStats = allRegion_mergedStats
## sensitivity model stats
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_mergedStats_DMP_analysis_dupCor_NeuN_sensitivity.rda')
sensitivityStats = allRegion_mergedStats[,25:44]
#
colnames(sensitivityStats) <- gsub("ALL_", "Sensitivity_", colnames(sensitivityStats) )
colnames(primaryStats) <- gsub("ALL_","Primary_", colnames(primaryStats) )

## Combine
allStats = cbind(primaryStats, sensitivityStats[rownames(primaryStats),] )
save(allStats, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda' )

### Cut FDR significant (<5%) main effect and interaction effect CpGs for supplemental tables
write.csv(allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05, ],row.names=F,file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/mainEffect_cpgs_fdr05_dupcor.csv')
write.csv(allStats[allStats$Primary_subset_interactionEffect_adj.P.Val<0.05,],row.names=F,file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/interaction_cpgs_fdr05_dupcor.csv')
##### Tidying up the stats object
tidyMergedStats= allStats %>%
  gather(key, value, 25:ncol(allStats))

###
tidyMergedStats$Model <- gsub("\\_[^\\_]*$","",tidyMergedStats[,'key'],)
tidyMergedStats$Statistic = gsub("(^.*_)","",tidyMergedStats[,'key'],)
###
tidyStats = tidyMergedStats[,-grep("key",colnames(tidyMergedStats))]  
tidyStats = spread(tidyStats,Statistic,value)
###
tidyStats$CellTypeAdjustment = jaffelab::ss(tidyStats$Model,"_",1)
tidyStats$Model = unlist(lapply(stringr::str_extract_all(tidyStats$Model, "[^_]+"), function(x) paste(x[-1],collapse ="_") ))
tidyStats$Interaction = jaffelab::ss(tidyStats$Model,"_",2)

allRegion_tidyStats = tidyStats
save(allRegion_tidyStats, file ='/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_allRegion.rda')
##
allRegion_DMP_MetaStats =group_by(allRegion_tidyStats, Interaction, CellTypeAdjustment, Model) %>% summarise( N_DMP_FDR_01 = sum(`adj.P.Val`<0.01),
												 N_DMP_FDR_05 = sum(`adj.P.Val`<0.05), 
												 DMP_Up_N =sum(t[adj.P.Val <0.05]>0), 
												 DMP_Down_N=sum(t[adj.P.Val <0.05]<0), 
												 DMP_Up_Percent = 100*DMP_Up_N/N_DMP_FDR_05, 
												 DMP_Median_EffectSize_Magnitude = median(abs(logFC[adj.P.Val <0.05])) ) %>% as.data.frame()
						
allRegion_DMP_MetaStats$DMP_Direction_P = NA
allRegion_DMP_MetaStats$DMP_Direction_P[!is.na(allRegion_DMP_MetaStats$DMP_Up_N)]  = apply(allRegion_DMP_MetaStats[!is.na(allRegion_DMP_MetaStats$DMP_Up_N),c('DMP_Up_N', 'DMP_Down_N')] ,1, function(x){ chisq.test(x)$p.value } )
write.csv(allRegion_DMP_MetaStats, '/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_caseControl_DMP_Results_Metastats_FDR05.csv',row.names=F)


## 
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
tidyMergedStats= regionSpecific_mergedStats %>%
  gather(key, value, 25:ncol(regionSpecific_mergedStats))

###
tidyMergedStats$Model <- gsub("\\_[^\\_]*$","",tidyMergedStats[,'key'],)
tidyMergedStats$Statistic = gsub("(^.*_)","",tidyMergedStats[,'key'],)
###
tidyStats = tidyMergedStats[,-grep("key",colnames(tidyMergedStats))]  
tidyStats = spread(tidyStats,Statistic,value)
###
tidyStats$Region = jaffelab::ss(tidyStats$Model,"_",1)
tidyStats$CellTypeAdjustment = jaffelab::ss(tidyStats$Model,"_",3)
tidyStats$Model = jaffelab::ss(tidyStats$Model,"_",2)
##
singleRegion_tidyStats = tidyStats
singleRegion_tidyStats[singleRegion_tidyStats$Region == "DLPFC" & singleRegion_tidyStats$Model == "subset" & singleRegion_tidyStats$CellTypeAdjustment =="NoAdj", 'adj.P.Val' ] <- p.adjust(singleRegion_tidyStats[singleRegion_tidyStats$Region == "DLPFC" & singleRegion_tidyStats$Model == "subset" & singleRegion_tidyStats$CellTypeAdjustment =="NoAdj", 'P.Value' ],method="fdr")
save(singleRegion_tidyStats, file ='/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_singleRegion.rda')

##
singleRegion_DMP_MetaStats =group_by(singleRegion_tidyStats, Region, Model, CellTypeAdjustment) %>% summarise( N_DMP_FDR_01 = sum(`adj.P.Val`<0.01),
												 N_DMP_FDR_05 = sum(`adj.P.Val`<0.05), 
												 N_DMP_FDR_10 = sum(`adj.P.Val`<0.10),
												 DMP_Up_N =sum(t[adj.P.Val <0.10]>0), 
												 DMP_Down_N=sum(t[adj.P.Val <0.10]<0), 
												 DMP_Up_Percent = 100*DMP_Up_N/N_DMP_FDR_10, 
												 DMP_Median_EffectSize_Magnitude = median(abs(logFC[adj.P.Val <0.10])) ) %>% as.data.frame()
						
singleRegion_DMP_MetaStats$DMP_Direction_P = NA
singleRegion_DMP_MetaStats$DMP_Direction_P[!is.na(singleRegion_DMP_MetaStats$DMP_Up_N) & (singleRegion_DMP_MetaStats$DMP_Up_N + singleRegion_DMP_MetaStats$DMP_Down_N >0 )]  = apply(singleRegion_DMP_MetaStats[!is.na(singleRegion_DMP_MetaStats$DMP_Up_N) & (singleRegion_DMP_MetaStats$DMP_Up_N + singleRegion_DMP_MetaStats$DMP_Down_N >0 ),c('DMP_Up_N','DMP_Down_N')] ,1, function(x){ chisq.test(x)$p.value } )

write.csv(singleRegion_DMP_MetaStats, '/dcl01/lieber/ajaffe/Steve/Alz/csvs/singleRegion_caseControl_DMP_Results_Metastats.csv',row.names=F)




