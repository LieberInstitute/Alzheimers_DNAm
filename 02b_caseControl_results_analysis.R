library(RColorBrewer)
library(pheatmap)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
write.csv(regionSpecific_mergedStats[rowSums(regionSpecific_mergedStats[,grep("subset_NoAdj_adj.P.Val", colnames(regionSpecific_mergedStats),value=T)]<0.05)>0,],file='singleRegion_caseControl_DMP_Results_FDR05_subset_notCellTypeAdjusted.csv')
###### Check correlation between different model pvalues
pval_modelCor = cor(-log10(regionSpecific_mergedStats[,grep("_P.Value", colnames(regionSpecific_mergedStats),value=TRUE )]) )
tstat_modelCor = cor(regionSpecific_mergedStats[,grep("_t", colnames(regionSpecific_mergedStats),value=TRUE )[-1] ]) 
lfc_modelCor = cor(regionSpecific_mergedStats[,grep("_logFC", colnames(regionSpecific_mergedStats),value=TRUE )]) 
####
pdf("plots/regionSpecific_model_clustering.pdf",h=10,w=10,onefile=TRUE)
pheatmap(pval_modelCor, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)
pheatmap(tstat_modelCor, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)
pheatmap(lfc_modelCor, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)	
pheatmap(pval_modelCor[grep("subset",rownames(pval_modelCor),value=TRUE ),grep("subset",rownames(pval_modelCor),value=TRUE )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)
pheatmap(tstat_modelCor[grep("subset",rownames(tstat_modelCor),value=TRUE ),grep("subset",rownames(tstat_modelCor),value=TRUE )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)
pheatmap(lfc_modelCor[grep("subset",rownames(lfc_modelCor),value=TRUE ),grep("subset",rownames(lfc_modelCor),value=TRUE )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)				
dev.off()
########### Get all stats
models_of_interest = c("Primary_subset_mainEffect_adj.P.Val", "Primary_subset_interactionEffect_adj.P.Val")
library(eulerr)
library(VennDiagram)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')
## sig
FDR05_main_vs_interaction =lapply( models_of_interest, function(model, thresh=0.05, col = 'Name') {
allStats[allStats[,model] < thresh,col] })

names(FDR05_main_vs_interaction) <- c("Cross-region", "Region-dependent")
fit = euler(FDR05_main_vs_interaction)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/main_effect_vs_interaction_effect_venn_FDR05.pdf',height=12,width=12)
plot(fit, fill_opacity = 0.3,quantities=T,cex=2 )
dev.off()

venn_fdr05 = venn.diagram(x = FDR05_main_vs_interaction,
							category.names = names(FDR05_main_vs_interaction),
							filename = NULL,
							fill = c('red', 'blue'),
							cat.just=list(c(0.9,1.5) , c(-0.8,5) ), cex=5, cat.cex=3
							)							
grid.draw(venn_fdr05)

	
pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/main_effect_vs_interaction_effect_vennDiagram_FDR05.pdf',height=12,width=12)
    grid.draw(venn_fdr05)
dev.off()
grep("DUSP", unique(allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05 & allStats$Primary_subset_interactionEffect_adj.P.Val<0.05,'UCSC_RefGene_Name']),value=T )

#####3
library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda',verbose=T)
###
library(ggplot2)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##Interaction
dat = regionSpecific_mergedStats[allStats[allStats$Primary_subset_interactionEffect_adj.P.Val<0.05,'Name'],c('CRB_subset_NoAdj_logFC','DLPFC_subset_NoAdj_logFC','HIPPO_subset_NoAdj_logFC','ERC_subset_NoAdj_logFC') ]
colnames(dat) <- jaffelab::ss(colnames(dat),"_",1)
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/GGally_pairwise_scatter_17056_interaction_cpgs_4brainRegions.pdf',height=10,width=10)
GGally::ggscatmat(dat, columns = 1:4 )
dev.off()

## Main
dat = regionSpecific_mergedStats[allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05,'Name'],c('CRB_subset_NoAdj_logFC','DLPFC_subset_NoAdj_logFC','HIPPO_subset_NoAdj_logFC','ERC_subset_NoAdj_logFC') ]
colnames(dat) <- jaffelab::ss(colnames(dat),"_",1)
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/GGally_pairwise_scatter_1494_main_cpgs_4brainRegions.pdf',height=10,width=10)
GGally::ggscatmat(dat, columns = 1:4 )
dev.off()

#tidyMergedStats = data.table::as.data.table(tidyMergedStats)
#tidyMergedStats[,c('Region','Model','CellTypeAdjustment','Statistic') := data.table::tstrsplit(tidyMergedStats[,'key'], "_", fixed=TRUE)]
#regionalTidyStats= tidyMergedStats

#tidyMergedStats$Region =  gsub("_.*","",tidyMergedStats[,'key'],)
#tidyMergedStats$Model =  jaffelab::ss(tidyMergedStats[,'key'], "_", 2)
#tidyMergedStats$CellTypeAdjustment =  jaffelab::ss(tidyMergedStats[,'key'], "_", 3)
#tidyMergedStats$Statistic =  gsub("\\_[^\\_]*$","",tidyMergedStats[1:10,'key'],)

#tidyMergedStats$statistic <- gsub("(^.*_)","",tidyMergedStats[1:10,'key'],)
#tidyMergedStats$model <- gsub("\\_[^\\_]*$","",tidyMergedStats[1:10,'key'],)

#tidyMergedStats2 = tidyMergedStats[,-grep("key",colnames(tidyMergedStats))]
#tidyMergedStats2 = spread(tidyMergedStats2,statistic,value)
#
#tidyStats = tidyMergedStats2
#tidyStats$Region = jaffelab::ss(tidyStats$model,"_",1)
#tidyStats$Model = unlist(lapply(stringr::str_extract_all(tidyStats$model, "[^_]+"), function(x) paste(x[-1],collapse ="_") ))
#save(tidyStats, file ='/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidy_merged_DMC_case_control_stats.rda')

#load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidy_merged_DMC_case_control_stats.rda')

####### Gene level venn diagrams
## Annotate each probe with nearest gene and distance
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
theTranscripts = annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene,codingOnly=TRUE)

an = bumphunter::annotateNearest(GRanges(seqnames =regionSpecific_mergedStats$chr_hg38, IRanges(start=regionSpecific_mergedStats$pos_hg38, end =regionSpecific_mergedStats$pos_hg38),strand="*"),  theTranscripts)

regionSpecific_mergedStats$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
regionSpecific_mergedStats$nearestGeneDist = an$dist

###
models_notAdj_subset = grep("subset", grep("NoAdj",models,value=T), value=T)
## sig
sig_p1e2_regional_noAdj_subset =lapply( models_notAdj_subset, function(model, thresh=1e-2, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e2_regional_noAdj_subset) <- c(jaffelab::ss(models_notAdj_subset,"_",1))

sig_p1e3_regional_noAdj_subset =lapply( models_notAdj_subset, function(model, thresh=1e-3, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e3_regional_noAdj_subset) <- c(jaffelab::ss(models_notAdj_subset,"_",1))

sig_p1e4_regional_noAdj_subset =lapply( models_notAdj_subset, function(model, thresh=1e-4, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e4_regional_noAdj_subset) <- c(jaffelab::ss(models_notAdj_subset,"_",1))

fit_2 <- euler(sig_p1e2_regional_noAdj_subset)
fit_3 <- euler(sig_p1e3_regional_noAdj_subset)
fit_4 <- euler(sig_p1e4_regional_noAdj_subset)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/regional_DMP_GeneLevel_Unadjusted_subset_vennDiagram_plots.pdf',height=12,width=12)
plot(fit_2, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e2_regional_noAdj_subset), order.by = "freq")

plot(fit_3, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e3_regional_noAdj_subset), order.by = "freq")

plot(fit_4, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e4_regional_noAdj_subset), order.by = "freq")

dev.off()  






### Adjustment#############
models_Adj_subset = grep("subset", grep("cellTypeAdj",models,value=T), value=T)
## Sig
sig_p1e2_regional_Adj_subset =lapply( models_Adj_subset, function(model, thresh=1e-2, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e2_regional_Adj_subset) <- c(jaffelab::ss(models_Adj_subset,"_",1))

sig_p1e3_regional_Adj_subset =lapply( models_Adj_subset, function(model, thresh=1e-3, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e3_regional_Adj_subset) <- c(jaffelab::ss(models_Adj_subset,"_",1))

sig_p1e4_regional_Adj_subset =lapply( models_Adj_subset, function(model, thresh=1e-4, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e4_regional_Adj_subset) <- c(jaffelab::ss(models_Adj_subset,"_",1))

fit_2 <- euler(sig_p1e2_regional_Adj_subset)
fit_3 <- euler(sig_p1e3_regional_Adj_subset)
fit_4 <- euler(sig_p1e4_regional_Adj_subset)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/regional_DMP_GeneLevel_CellTypeAdjusted_subset_vennDiagram_plots.pdf',height=12,width=12)
plot(fit_2, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e2_regional_Adj_subset), order.by = "freq")

plot(fit_3, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e3_regional_Adj_subset), order.by = "freq")

plot(fit_4, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple') )
upset(fromList(sig_p1e4_regional_Adj_subset), order.by = "freq")

dev.off() 



#######
library(VennDiagram)
subset_models=grep("subset_a", grep("_adj.P.Val", colnames(mergedStats),value=TRUE ), value=TRUE)
full_models=grep("full_a", grep("_adj.P.Val", colnames(mergedStats),value=TRUE ), value=TRUE)
subset05=lapply( subset_models, function(model) {
mergedStats[mergedStats[,model] <0.05,'Name']
})
subset10=lapply( subset_models, function(model) {
mergedStats[mergedStats[,model] <0.10,'Name']
})

full05=lapply( full_models, function(model) {
mergedStats[mergedStats[,model] <0.05,'Name']
})
full10=lapply( full_models, function(model) {
mergedStats[mergedStats[,model] <0.10,'Name']
})

names(subset05) <- jaffelab::ss(subset_models,"_",1)
names(subset10) <- jaffelab::ss(subset_models,"_",1)
names(full05) <- jaffelab::ss(full_models,"_",1)
names(full10) <- jaffelab::ss(full_models,"_",1)

subset_venn05 = venn.diagram(x = subset05,
							category.names = names(subset05),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow')
							)
subset_venn10 = venn.diagram(x = subset10,
							category.names = names(subset05),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow')
							)
full_venn05 = venn.diagram(x = full05,
							category.names = names(full05),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow')
							)
full_venn10 = venn.diagram(x = full10,
							category.names = names(full05),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow')
							)							
pdf(file='plots/Sig_DMC_regionalSubsetModels_vennDiagram.pdf')
    grid.draw(subset_venn05)
	grid.newpage()
    grid.draw(subset_venn10)
dev.off()
pdf(file='plots/Sig_DMC_regionalFullModels_vennDiagram.pdf')
    grid.draw(full_venn05)
	grid.newpage()
    grid.draw(full_venn10)
dev.off()
