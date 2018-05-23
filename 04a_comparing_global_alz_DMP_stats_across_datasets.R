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
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')

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
ii = rownames(allRegion_mergedStats)
multisetStats = cbind(allRegion_mergedStats,
					  Watson_mergedStats[match(ii, rownames(Watson_mergedStats)),],
					  Lunnon_mergedStats[match(ii, rownames(Lunnon_mergedStats)),],
					  Smith_mergedStats[match(ii, rownames(Smith_mergedStats)),])
### Using LFC
datLFC = multisetStats[,c('Name',
					   "ALL_subset_mainEffect_logFC",
					    "ALL_subset_mainEffect_adj.P.Val",
					   'STG_Dx_Unadjusted_logFC_Watson2016',
					   'ALL_main_Dx_logFC_Lunnon2014',
					   'STG_Dx_logFC_Lunnon2014',
					   'CRB_Dx_logFC_Lunnon2014',
					   'ERC_Dx_logFC_Lunnon2014',
					   'PFC_Dx_logFC_Lunnon2014')]
cor(datLFC[datLFC$ALL_subset_mainEffect_adj.P.Val<0.05,grep("logFC",colnames(datLFC))],use='pairwise.complete.obs')					   					   
colnames(datLFC)[4:ncol(datLFC)] = paste0(jaffelab::ss(colnames(datLFC)[4:ncol(datLFC)],"_",1),"_",gsub(".*_","",colnames(datLFC)[4:ncol(datLFC)]))			
					   
dat=tidyr::gather(datLFC,PublicSet,logFC,4:ncol(datLFC) )
pdf('plots/mainEffect_cpg_effectSizes_vs_public_data_effectSizes.pdf',height=8,width=14)	 
ggplot(subset(dat,ALL_subset_mainEffect_adj.P.Val<0.05 ), aes(x=ALL_subset_mainEffect_logFC, y=logFC))+facet_wrap(~PublicSet,nrow=1) +geom_point() + geom_smooth(se=F,method='lm') + labs(x='Our Data Effect Size', y='Public Data Effect Size')
dev.off()
### Using Pvalue

lun = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC) )
fisher.test(lun)

smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,
multisetStats$ALL_main_Braak_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Braak_logFC_Smith2018) == sign(multisetStats$ALL_subset_mainEffect_logFC) )
fisher.test(smith)

smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Dx_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) )

smith = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Braak_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Braak_logFC_Smith2018) == sign(multisetStats$ALL_subset_mainEffect_logFC) )


watson = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$STG_Dx_Unadjusted_P.Value_Watson2016<0.05 & sign(multisetStats$STG_Dx_Unadjusted_logFC_Watson2016) == sign(multisetStats$ALL_subset_mainEffect_logFC) )

###
anySet = watson = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05, 
(multisetStats$STG_Dx_Unadjusted_P.Value_Watson2016<0.05 & sign(multisetStats$STG_Dx_Unadjusted_logFC_Watson2016) == sign(multisetStats$ALL_subset_mainEffect_logFC) ) | 
(multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC)) | 
(multisetStats$ALL_main_Dx_P.Value_Smith2018<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) ) ) 

topHits = multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & 
(multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC))  &multisetStats$ALL_main_Braak_P.Value_Smith2018<0.05& sign(multisetStats$ALL_main_Dx_logFC_Smith2018) != sign(multisetStats$ALL_subset_mainEffect_logFC) 
multisetStats[topHits,'within10kb_geneSymbol_gencode_hg38']
####### Saving out replication results

write.csv(multisetStats[multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC),],file='csvs/SupplementalTable_dmp_replicated_nomSig_dirConsistent_in_lunnon2014.csv',row.names=F )

###
j = table(multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05,multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025  )

write.csv(multisetStats[multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05&multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025,], file='csvs/dmp_replicated_nomSig_interaction_in_lunnon2014.csv',row.names=F )
