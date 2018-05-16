library(RColorBrewer)
library(pheatmap)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(eulerr)
library(UpSetR)

col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
#regionSpecific_mergedStats$DLPFC_subset_NoAdj_adj.P.Val <- p.adjust(regionSpecific_mergedStats$DLPFC_subset_NoAdj_P.Value,method='fdr')

#### Some Venn Diagrams 

models=grep("P.Value", colnames(regionSpecific_mergedStats),value=TRUE )
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

### FDR Significant Venn Diagrams
models_notAdj_subset_FDR = gsub("P.Value", "adj.P.Val", models_notAdj_subset)

sig_FDR10_regional_noAdj_subset =lapply( models_notAdj_subset_FDR, function(model, thresh=.10, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_FDR10_regional_noAdj_subset) <- c(jaffelab::ss(models_notAdj_subset_FDR,"_",1))
sig_FDR10_regional_noAdj_subset = sig_FDR10_regional_noAdj_subset[lengths(sig_FDR10_regional_noAdj_subset)!=0]


fit_2 <- euler(sig_p1e2_regional_noAdj_subset)
fit_3 <- euler(sig_p1e3_regional_noAdj_subset)
fit_4 <- euler(sig_p1e4_regional_noAdj_subset)
fit_FDR10 <- euler(sig_FDR10_regional_noAdj_subset)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/regional_DMP_Unadjusted_subset_vennDiagram_plots.pdf',height=12,width=12)
plot(fit_2, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'),cex=1.4)
upset(fromList(sig_p1e2_regional_noAdj_subset), order.by = "freq", text.scale=3 )

plot(fit_3, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'), cex=1.4 )
upset(fromList(sig_p1e3_regional_noAdj_subset), order.by = "freq",text.scale=3)

plot(fit_4, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'), cex=1.4 )
upset(fromList(sig_p1e4_regional_noAdj_subset), order.by = "freq",text.scale=3)

plot(fit_FDR10, fill_opacity = 0.3,quantities=T, fill=c('red','blue'), cex=1.4 )

dev.off()  

### Adjustment
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

### FDR Significant Venn Diagrams
models_Adj_subset_FDR = gsub("P.Value", "adj.P.Val", models_Adj_subset)

sig_FDR10_regional_Adj_subset =lapply( models_Adj_subset_FDR, function(model, thresh=.10, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_FDR10_regional_Adj_subset) <- c(jaffelab::ss(models_Adj_subset_FDR,"_",1))
sig_FDR10_regional_Adj_subset = sig_FDR10_regional_Adj_subset[lengths(sig_FDR10_regional_Adj_subset)!=0]

fit_2 <- euler(sig_p1e2_regional_Adj_subset)
fit_3 <- euler(sig_p1e3_regional_Adj_subset)
fit_4 <- euler(sig_p1e4_regional_Adj_subset)
fit_FDR10 <- euler(sig_FDR10_regional_Adj_subset)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/regional_DMP_CellTypeAdjusted_subset_vennDiagram_plots.pdf',height=12,width=12)
plot(fit_2, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'),cex=1.4 )
upset(fromList(sig_p1e2_regional_Adj_subset), order.by = "freq",text.scale=2)

plot(fit_3, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'),cex=1.4 )
upset(fromList(sig_p1e3_regional_Adj_subset), order.by = "freq",text.scale=2)

plot(fit_4, fill_opacity = 0.3,quantities=T, fill=c('red','blue','green','purple'),cex=1.4 )
upset(fromList(sig_p1e4_regional_Adj_subset), order.by = "freq",text.scale=2)

plot(fit_FDR10, fill_opacity = 0.3,quantities=T, fill=c('red','blue'),cex=1.4 )

dev.off() 



## Sig
sig_p1e2_regional_Adj_subset =lapply( models_Adj_subset, function(model, thresh=1e-2, col = 'Name') {
regionSpecific_mergedStats[regionSpecific_mergedStats[,model] < thresh,col] })
names(sig_p1e2_regional_Adj_subset) <- c(jaffelab::ss(models_Adj_subset,"_",1))

### allRegion venn diagrams 
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')

models=grep("P.Value", colnames(allStats),value=TRUE )

###
models_subset = grep("subset", models, value=T)
## sig
sig_p1e2_allRegion_subset =lapply( models_subset, function(model, thresh=1e-2, col = 'Name') {
allStats[allStats[,model] < thresh,col] })
names(sig_p1e2_allRegion_subset) <- c(paste0(jaffelab::ss(models_subset,"_",1),"_",jaffelab::ss(models_subset,"_",3)) )

sig_p1e3_allRegion_subset =lapply( models_subset, function(model, thresh=1e-3, col = 'Name') {
allStats[allStats[,model] < thresh,col] })
names(sig_p1e3_allRegion_subset) <- c(paste0(jaffelab::ss(models_subset,"_",1),"_",jaffelab::ss(models_subset,"_",3)) )

sig_p1e4_allRegion_subset =lapply( models_subset, function(model, thresh=1e-4, col = 'Name') {
allStats[allStats[,model] < thresh,col] })
names(sig_p1e4_allRegion_subset) <- c(paste0(jaffelab::ss(models_subset,"_",1),"_",jaffelab::ss(models_subset,"_",3)) )

fit_2 <- euler(sig_p1e2_allRegion_subset)
fit_3 <- euler(sig_p1e3_allRegion_subset)
fit_4 <- euler(sig_p1e4_allRegion_subset)

pdf(file='/dcl01/lieber/ajaffe/Steve/Alz/plots/allRegion_DMP_subset_vennDiagram_plots.pdf',height=12,width=12)
plot(fit_2, fill_opacity = 0.3,quantities=T, fill=c('#ca5d57','#399283','#9f6cbd','#989c3f'),cex=1.4, main = "P<1e-2" )
upset(fromList(sig_p1e2_allRegion_subset), order.by = "freq", text.scale = 2)

plot(fit_3, fill_opacity = 0.3,quantities=T, fill=c('#ca5d57','#399283','#9f6cbd','#989c3f'),cex=1.4, main = "P<1e-3" )
upset(fromList(sig_p1e3_allRegion_subset), order.by = "freq", text.scale = 2)

plot(fit_4, fill_opacity = 0.3,quantities=T, fill=c('#ca5d57','#399283','#9f6cbd','#989c3f'),main = "P<1e-4",cex=1.4 )
upset(fromList(sig_p1e4_allRegion_subset), order.by = "freq", text.scale = 2)

dev.off() 

### overlap between 'main effect' cpgs and interaction cpgs

