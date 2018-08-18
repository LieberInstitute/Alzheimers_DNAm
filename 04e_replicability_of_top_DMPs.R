library(dplyr)
library(tidyr)
check_replication= function(discovery_names = a1,
							discovery_effect_size= a2,
							replication_names = a3,
							replication_p = a4,
							replication_effect_size = a5) {
shared_names = unique(intersect(discovery_names, replication_names))
disc_i=match(shared_names, discovery_names)
rep_i=match(shared_names, replication_names)
							
N = sum(replication_p[rep_i] <0.05 & sign(discovery_effect_size[disc_i]) == sign(replication_effect_size[rep_i]),na.rm=F )

percent = N*100/length(shared_names)
return( c(N,percent) )
}



###Check replicability of our top genes 

## Load our statistics
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_singleRegion.rda')

## Load watson et al 2016: STG
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE76105/GSE76105_Watson_mergedStats.rda')
colnames(Watson_mergedStats) = paste0(colnames(Watson_mergedStats), ": ", "Watson2016")

## Load lunnon et al 2014: Multiregion
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_Lunnon_mergedStats.rda')
colnames(Lunnon_mergedStats) = paste0(colnames(Lunnon_mergedStats), ": ", "Lunnon2014")

##
sigSingleRegion_tidyStats = singleRegion_tidyStats[ singleRegion_tidyStats[,'adj.P.Val']<0.1, ]
										 
singleRegion_DMP_Replication=group_by(sigSingleRegion_tidyStats,Region, Model, CellTypeAdjustment) %>% summarise(
FDR_10=n() , 
												 Lunnon_PFC_Dx_N= check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`PFC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`PFC_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_PFC_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`PFC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`PFC_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_STG_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`STG_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`STG_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_STG_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`STG_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`STG_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_ERC_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ERC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ERC_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_ERC_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ERC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ERC_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_CRB_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`CRB_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`CRB_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_CRB_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`CRB_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`CRB_Dx_logFC: Lunnon2014`)[2],
												 
												 Watson_STG_Dx_unAdj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_unAdj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[2],
												 
												 Watson_STG_Dx_Adj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_Adj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[2]) %>% as.data.frame()
												 

top100SingleRegion_tidyStats = group_by(singleRegion_tidyStats, Region, Model, CellTypeAdjustment) %>% top_n(n = 100, wt = `P.Value`)
												 
singleRegion_DMP_Replication_top100=group_by(top100SingleRegion_tidyStats,Region, Model, CellTypeAdjustment) %>% summarise(
N_CpG=n() , 
												 Lunnon_PFC_Dx_N= check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`PFC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`PFC_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_PFC_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`PFC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`PFC_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_STG_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`STG_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`STG_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_STG_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`STG_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`STG_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_ERC_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ERC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ERC_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_ERC_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ERC_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ERC_Dx_logFC: Lunnon2014`)[2],
												 
												 Lunnon_CRB_Dx_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`CRB_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`CRB_Dx_logFC: Lunnon2014`)[1],
												 
												 Lunnon_CRB_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`CRB_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`CRB_Dx_logFC: Lunnon2014`)[2],
												 
												 Watson_STG_Dx_unAdj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_unAdj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[2],
												 
												 Watson_STG_Dx_Adj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_Adj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[2]) %>% as.data.frame()

singleRegion_DMP_Replication_top100 = singleRegion_DMP_Replication_top100[ ,-grep("_N", colnames(singleRegion_DMP_Replication_top100))]

openxlsx::write.xlsx(list(FDR10_replication=singleRegion_DMP_Replication,top100_replication=singleRegion_DMP_Replication_top100), file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/LIBD_Discovery_Replication_Results.xlsx')		


## Load our statistics: ALL REGION
			
## Load our statistics
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_allRegion.rda')

##
sigAllRegion_tidyStats = allRegion_tidyStats[ allRegion_tidyStats[,'adj.P.Val']<0.1, ]
										 
allRegion_DMP_Replication=sigAllRegion_tidyStats%>%filter(Interaction=="mainEffect")%>%group_by(Sensitivity, Model) %>% summarise(
FDR_10=n() , 
												 Lunnon_Main_Dx_N= check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ALL_main_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ALL_main_logFC: Lunnon2014`)[1],
												 
												 Lunnon_Main_Dx_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Lunnon_mergedStats$`ALL_main_Dx_P.Value: Lunnon2014`,
												 replication_effect_size=Lunnon_mergedStats$`ALL_main_Dx_logFC: Lunnon2014`)[2],
												 
												 Watson_STG_Dx_unAdj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_unAdj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_Unadjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_Unadjusted_logFC: Watson2016`)[2],
												 
												 Watson_STG_Dx_Adj_N = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[1],
												 
												 Watson_STG_Dx_Adj_Percent = check_replication( discovery_names=Name,
												 discovery_effect_size=logFC,
												 replication_names=Name,
												 replication_p=Watson_mergedStats$`STG_Dx_cellCompAdjusted_P.Value: Watson2016`,
												 replication_effect_size=Watson_mergedStats$`STG_Dx_cellCompAdjusted_logFC: Watson2016`)[2]) %>% as.data.frame()			