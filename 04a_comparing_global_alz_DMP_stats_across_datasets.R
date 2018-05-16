## Comparing global alzheimer statistics 
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
## Load our statistics
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_singleRegion.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')

## Load watson et al 2016: STG
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE76105/GSE76105_Watson_mergedStats.rda')
colnames(Watson_mergedStats) = paste0(colnames(Watson_mergedStats), "_", "Watson2016")

## Load lunnon et al 2014: Multiregion
load('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_Lunnon_mergedStats.rda')
colnames(Lunnon_mergedStats) = paste0(colnames(Lunnon_mergedStats), "_", "Lunnon2014")

##
ii = rownames(allRegion_mergedStats)
multisetStats = cbind(allRegion_mergedStats,
					  Watson_mergedStats[match(ii, rownames(Watson_mergedStats)),],
					  Lunnon_mergedStats[match(ii, rownames(Lunnon_mergedStats)),])
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
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/mainEffect_cpg_effectSizes_vs_public_data_effectSizes.pdf',height=8,width=14)	 
ggplot(subset(dat,ALL_subset_mainEffect_adj.P.Val<0.05 ), aes(x=ALL_subset_mainEffect_logFC, y=logFC))+facet_wrap(~PublicSet,nrow=1) +geom_point() + geom_smooth(se=F,method='lm') + labs(x='Our Data Effect Size', y='Public Data Effect Size')
dev.off()
### Using Pvalue

k = table(multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05,multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC) )
fisher.test(k)
write.csv(multisetStats[multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC),],file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/dmp_replicated_nomSig_dirConsistent_in_lunnon2014.csv',row.names=F )
###
j = table(multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05,multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025  )
write.csv(multisetStats[multisetStats$ALL_subset_interactionEffect_adj.P.Val<0.05&multisetStats$ALL_int_Dx_P.Value_Lunnon2014<0.025,], file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/dmp_replicated_nomSig_interaction_in_lunnon2014.csv',row.names=F )


##
sig=multisetStats[multisetStats$ALL_subset_mainEffect_adj.P.Val<0.05 & multisetStats$ALL_main_Dx_P.Value_Lunnon2014<0.05 & sign(multisetStats$ALL_main_Dx_logFC_Lunnon2014) == sign(multisetStats$ALL_subset_mainEffect_logFC),"Name"]
allCpG=unique(multisetStats$Name)
sigCpg= sig
gstGO <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=allCpG, collection="GO")
gstGO[grep("cell adhesion",gstGO$Term),]
gstGO[order(gstGO$P.DE),][1:10,]

gstKEGG <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=allCpG, collection="KEGG")
gstKEGG[order(gstKEGG$P.DE),][1:10,]

#####

datPvalue = multisetStats[,c('Name',
					   "ALL_subset_mainEffect_P.Value",
					   "ALL_subset_mainEffect_adj.P.Val",
					   "ALL_subset_interactionEffect_P.Value",
					   "ALL_subset_interactionEffect_adj.P.Val",				   
					   'STG_Dx_Unadjusted_P.Value_Watson2016',
					   'ALL_main_Dx_P.Value_Lunnon2014',
					   'STG_Dx_P.Value_Lunnon2014',
					   'CRB_Dx_P.Value_Lunnon2014',
					   'ERC_Dx_P.Value_Lunnon2014',
					   'PFC_Dx_P.Value_Lunnon2014')]

tab = table(datPvalue$ALL_subset_mainEffect_adj.P.Val<0.05,matrixStats::rowMins(as.matrix(datPvalue[,c('STG_Dx_Unadjusted_P.Value_Watson2016',
					   'ALL_main_Dx_P.Value_Lunnon2014',
					   'STG_Dx_P.Value_Lunnon2014',
					   'CRB_Dx_P.Value_Lunnon2014',
					   'ERC_Dx_P.Value_Lunnon2014',
					   'PFC_Dx_P.Value_Lunnon2014')]) )<1e-2 )
tab
fisher.test(tab)

tab = table(datPvalue$ALL_subset_interactionEffect_adj.P.Val<0.05,matrixStats::rowMins(as.matrix(datPvalue[,c('STG_Dx_Unadjusted_P.Value_Watson2016',
					   'ALL_main_Dx_P.Value_Lunnon2014',
					   'STG_Dx_P.Value_Lunnon2014',
					   'CRB_Dx_P.Value_Lunnon2014',
					   'ERC_Dx_P.Value_Lunnon2014',
					   'PFC_Dx_P.Value_Lunnon2014')]) )<0.05 )
tab
fisher.test(tab)


##					  
library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')

adj_pvals = as.matrix(multisetStats[,grep("adj.P",colnames(multisetStats),value=T)]
adj_pvals = adj_pvals[,-grep("int",colnames(adj_pvals) )]
sigCpgs = which(rowSums(adj_pvals<0.1)>0)

pvals = -log10(as.matrix(multisetStats[,grep("P.Value",colnames(multisetStats),value=T)] ))
drop_full_pvalue_cols = grep("full",colnames(pvals))

effect_size = as.matrix(multisetStats[,grep("logFC",colnames(multisetStats),value=T)] )
drop_full_effect_cols = grep("full",colnames(effect_size))
pdf("/dcl01/lieber/ajaffe/Steve/Alz/plots/pheatmap_lunnon_watson_lieber_pvalue_gl_correlation.pdf",h=30,w=30,onefile=TRUE)
pheatmap(cor(pvals[,-drop_full_pvalue_cols], use='pairwise.complete.obs'), 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=30)
pheatmap(effect_size[sigCpgs,-drop_full_effect_cols], 
		cluster_rows=T, 
		cluster_cols=T,
		fontsize=30)			
dev.off()
	
#####
library(bacon)
library(BiocParallel)
register(MulticoreParam(1))

multisetStats = multisetStats[!(is.na(multisetStats[,'STG_Dx_Unadjusted_adj.P.Val: Watson2016'])),]

tstat=multisetStats[,gsub("SE","t",grep("SE",colnames(multisetStats),value=T))]


effect_size=multisetStats[,gsub("SE","logFC",grep("SE",colnames(multisetStats),value=T))]
se=multisetStats[,grep("SE",colnames(multisetStats),value=T)]
cols <- gsub("_t.*","",colnames(tstat))
colnames(tstat) = colnames(effect_size) = colnames(se) = cols

library(metafor)
res.fe <- metafor::rma(yi=t(effect_size[1, ]), sei=t(se[1, ]), method="FE")
forest.rma(res.fe, slab=gsub("_eset$","",colnames(effect_size) ))

#exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
#
#Make plot
#plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1)
#abline(0,1)


bakin = bacon( NULL,effect_size[1:2000,], se[1:2000,] )

bakin_meta <- meta(bakin)

es <- replicate(6, rnormmix(2000, c(0.9, 0, 1, 0, 4, 1)))
se1 <- replicate(6, 0.8*sqrt(4/rchisq(2000,df=4)))
bc <- bacon(NULL, es, se1)
mbc <- meta(bc)
BiocParallel::register(BiocParallel::MulticoreParam(1))

metafor::rma()
#The resulting coefficient values (interpreted as log hazard ratios) and standard errors were combined using the R software package metafor [34]. The meta-analysis was carried out with the R command rma (with arguments method=“FE” to get fixed effects estimates). The forest plots were created using the R function forest (with argument atransf=exp to exponentiate the estimate of the log hazard ratios).	
