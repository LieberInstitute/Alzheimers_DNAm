library(ggplot2)
library(scales)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5) ))
## load data
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda',verbose=T)
primaryStats = allRegion_mergedStats

## Test for enrichment in CpG Islands
allRegion_mergedStats$Relation_to_Island = plyr::revalue(allRegion_mergedStats$Relation_to_Island, c("N_Shore"="Shore", "S_Shore"="Shore", "N_Shelf"="Shelf","S_Shelf"="Shelf") )
chisq.test ( table(allRegion_mergedStats$Relation_to_Island, allRegion_mergedStats$ALL_subset_mainEffect_logFC < 0 & allRegion_mergedStats$ALL_subset_mainEffect_adj.P.Val <0.05) )
chisq.test (table(allRegion_mergedStats$Relation_to_Island, allRegion_mergedStats$ALL_subset_mainEffect_logFC > 0 & allRegion_mergedStats$ALL_subset_mainEffect_adj.P.Val <0.05) )
chisq.test (table(allRegion_mergedStats$Relation_to_Island, allRegion_mergedStats$ALL_subset_mainEffect_adj.P.Val <0.05) )

#Make plot
load('rdas/tidyStats_caseControl_DMC_allRegion.rda')
dat_sig = allRegion_tidyStats[allRegion_tidyStats$Interaction=="mainEffect" & allRegion_tidyStats$CellTypeAdjustment=="Primary" &allRegion_tidyStats$Model=="subset_mainEffect", ]
dat_sig$Relation_to_Island = plyr::revalue(dat_sig$Relation_to_Island, c("N_Shore"="Shore", "S_Shore"="Shore", "N_Shelf"="Shelf","S_Shelf"="Shelf") )

dat1 = data.frame(X="Background", Y = dat_sig$Relation_to_Island )
dat2 = data.frame(X="All Sig", Y = dat_sig[dat_sig$adj.P.Val<0.05,'Relation_to_Island'] )
dat3 = data.frame(X="Hyper Sig", Y = dat_sig[dat_sig$adj.P.Val<0.05 & dat_sig$logFC>0,'Relation_to_Island'] )
dat4 = data.frame(X="Hypo Sig", Y = dat_sig[dat_sig$adj.P.Val<0.05& dat_sig$logFC<0,'Relation_to_Island'] )
dat= rbind(dat1,dat2,dat3,dat4) 
library(dplyr)
dat = group_by(dat, X,Y) %>% summarise(count=n() )
dat$Y=factor(dat$Y,levels=c('OpenSea','Shelf','Shore','Island') )

cpg_main= ggplot(dat,aes(x = X, y =count,fill = Y)) + 
    geom_bar(position = "fill",stat = "identity") + labs(x="",y="",fill="") +
     geom_bar(position = position_fill(), stat = "identity") +scale_fill_grey() + scale_y_continuous(labels = percent_format()  )  + theme_bw(base_size=40) + theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(cpg_main, file='plots/SupplementalFigure_maineffect_dmp_enrichment_in_cpg.pdf',height=14,width=10,useDingbats=FALSE )

## sensitivity model stats
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor_NeuN_sensitivity.rda')
sensitivityStats = allRegion_mergedStats[,40:59]
### Scatterplot of effect sizes?
cor.test(primaryStats[,'ALL_subset_mainEffect_t'], sensitivityStats[rownames(primaryStats),"ALL_subset_mainEffect_t"])

pdf('plots/SupplementalFigure_main_effect_subsetModel_unadjusted_vs_adjusted_cell_type_senstivity.pdf',useDingbats=FALSE)

par(mar=c(5,6,2,2))
plot(primaryStats[,'ALL_subset_mainEffect_t'], 
	 sensitivityStats[rownames(primaryStats),"ALL_subset_mainEffect_t"], pch = 21, bg="grey",
	xlab="T-statistic without cell type adjustment", ylab="T-statistic with cell type adjustment",
	xlim=c(-8,8), ylim = c(-8,8),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-8,8), y = c(-8,8),col="red")
dev.off()

table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 sensitivity=sensitivityStats[rownames(primaryStats),"ALL_subset_mainEffect_adj.P.Val"]<0.05)
	 
### ordinal model
cor.test(primaryStats[,'ALL_subset_mainEffect_t'], primaryStats[rownames(primaryStats),"ALL_ordinal_mainEffect_t"])

pdf('plots/SupplementalFigure_main_effect_subsetModel_primary_vs_ordinal_secondary_model_tstatistics.pdf',useDingbats=FALSE)
par(mar=c(5,6,2,2))
plot(primaryStats[,'ALL_subset_mainEffect_t'], 
	 primaryStats[rownames(primaryStats),"ALL_ordinal_mainEffect_t"], pch = 21, bg="grey",
	xlab="Primary model t-statistic", ylab="Ordinal model t-statistic",
	xlim=c(-8,8), ylim = c(-8,8),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-8,8), y = c(-8,8),col="red")
dev.off()

table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 ordinal=primaryStats[rownames(primaryStats),"ALL_ordinal_mainEffect_P.Value"]<0.01)

	 ### APOE Sensitivity
load('rdas/APOE4_sensitivity_allRegion_mergedStats_DMP_analysis_dupCor.rda',verbose=T)

cor.test(primaryStats[,'ALL_subset_mainEffect_t'], APOE_allRegionStats[rownames(primaryStats),"ALL_subset_mainEffect_t"])

png('plots/SupplementalFigure_main_effect_subsetModel_primary_vs_APOE_Sensitivitiy_model_tstatistics.png')
par(mar=c(5,6,2,2))
plot(primaryStats[,'ALL_subset_mainEffect_t'], 
	 APOE_allRegionStats[rownames(primaryStats),"ALL_subset_mainEffect_t"],
	 , pch = 21, bg="grey",
	xlab="Primary model t-statistic", ylab="Model with APOE4 dosage t-statistic",
	xlim=c(-8,8), ylim = c(-8,8),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-8,8), y = c(-8,8),col="red")
dev.off()

table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 ordinal=APOE_allRegionStats[rownames(primaryStats),"ALL_subset_mainEffect_P.Value"]<0.01)


fisher.test(table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 ordinal=APOE_allRegionStats[rownames(primaryStats),"ALL_subset_mainEffect_adj.P.Val"]<0.05))

############ Posthoc effects size plots
library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda',verbose=T)
load('rdas/merged_DMP_regionSpecific_caseControl_stats.rda',verbose=T)
###
library(ggplot2)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##Interaction
dat = regionSpecific_mergedStats[allRegion_mergedStats[allRegion_mergedStats$ALL_subset_interactionEffect_adj.P.Val<0.05,'Name'],c('CRB_subset_NoAdj_logFC','DLPFC_subset_NoAdj_logFC','HIPPO_subset_NoAdj_logFC','ERC_subset_NoAdj_logFC') ]
colnames(dat) <- jaffelab::ss(colnames(dat),"_",1)
pdf('plots/SupplementalFigure_GGally_pairwise_scatter_regionDependent_DMPs_4brainRegions.pdf',height=10,width=10,useDingbats=FALSE)
GGally::ggscatmat(dat, columns = 1:4 )
dev.off()

## Main
dat = regionSpecific_mergedStats[allRegion_mergedStats[allRegion_mergedStats$Primary_subset_mainEffect_adj.P.Val<0.05,'Name'],c('CRB_subset_NoAdj_logFC','DLPFC_subset_NoAdj_logFC','HIPPO_subset_NoAdj_logFC','ERC_subset_NoAdj_logFC') ]
colnames(dat) <- jaffelab::ss(colnames(dat),"_",1)
pdf('plots/SupplementalFigure_GGally_pairwise_scatter_crossRegion_DMPs_4brainRegions.pdf',height=10,width=10,useDingbats=FALSE)
GGally::ggscatmat(dat, columns = 1:4 )
dev.off()

##############	 
	 
#Histogram for effect size 

######### 5% FDR Gene
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda',verbose=T)

histo_gene <- ggplot() + geom_histogram(data = allRegion_mergedStats[which(allRegion_mergedStats$ALL_subset_mainEffect_logFC < 0 & allRegion_mergedStats$ALL_subset_mainEffect_adj.P.Val <0.05), ],
                          aes(x = ALL_subset_mainEffect_logFC, y = ..count..),
						  binwidth=.01,
						  fill = "blue",
						  colour = "black") +
           geom_histogram(data = allRegion_mergedStats[which(allRegion_mergedStats$ALL_subset_mainEffect_logFC > 0 & allRegion_mergedStats$ALL_subset_mainEffect_adj.P.Val <0.05), ],
                          aes(x = ALL_subset_mainEffect_logFC, y = ..count..),
						  binwidth=.01,
						  fill = "red",
						  colour = "black") + 
			labs(x = "Change in DNAm",
				 y = "Number of Sites")	+
			theme_bw(base_size = 45) + 
		   theme(axis.title=element_text(size=50),
				 plot.title = element_text(hjust = 0.5),
				 panel.grid.major = element_blank(), 
				 panel.grid.minor = element_blank()) +
scale_x_continuous(limits = c(-0.2, 0.2))  + scale_y_continuous(limits = c(0, 150), expand = c(0, 0)) + 
 annotate("text", x = -0.10, y = 130, label = paste0("367 \nhypomethylated"),size=8) +
 annotate("text", x = 0.135, y = 130, label = paste0("491 \nhypermethylated"),size=8)
 
ggsave(histo_gene, file = 'plots/Figure_Histogram_of_DNAm_Percent_Change_MainSubsetDupCorModel.pdf',height=10,width=10,useDingbats=FALSE)
		   
##################### Aging analysis #####################
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/aging_controlOnly_mergedStats.rda')

##
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')
load('rdas/merged_DMP_regionSpecific_caseControl_stats.rda')

### Scatterplot of effect sizes?
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_t"])

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_mainEffect_P.Value']),-log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]) )

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_mainEffect_P.Value']),-log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]) )



cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"DLPFC_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"HIPPO_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ERC_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"CRB_Control_Aging_t"])

###
cor.test(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]) )

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]) )

#
png('plots/SupplementalFigure_aging_controlsOnly_vs_AD_Effect.png',res=300,units='in',height=7,width=7)
par(mar=c(5,6,2,2))
plot(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], 
	 Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_t"], pch = 21, bg="grey",
	xlab="Cross-region Alzheimer's T-statistic", ylab="Aging T-statistic (Controls Only)",
	xlim=c(-12,12), ylim = c(-12,12),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-12,12), y = c(-12,12),col="red")
dev.off()

png('plots/SupplementalFigure_aging_controlsOnly_vs_AD_Effect_regionDependentModel.png',res=300,units='in',height=7,width=7)
par(mar=c(5,6,2,2))
plot(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), 
	 -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]), pch = 21, bg="grey",
	xlab="-log10(Alzheimer's Interaction Pvalue)", ylab="-log10(Aging Interaction Pvalue)",
	xlim=c(0,20), ylim = c(0,20),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(0,20), y = c(0,20),col="red")
dev.off()

table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 sensitivity=sensitivityStats[rownames(primaryStats),"ALL_subset_mainEffect_adj.P.Val"]<0.05)
