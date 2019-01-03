library(ggplot2)
library(ggrepel)
library(jaffelab)
library(minfi)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

### load pca (n380)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/pca_alz.rda')


## load processed data
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)

###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  


#################### Checking correlation between PC1 and estimated NeuN+ prop. ####################
pd$PC1 = pca[pd$Sample_Name,1]
cor.test(pd$NeuN_pos, pd$PC1)

cor.test(pd[pd$Region=="CRB",'NeuN_pos'], pd[pd$Region=="CRB",'PC1'])
cor.test(pd[pd$Region=="DLPFC",'NeuN_pos'], pd[pd$Region=="DLPFC",'PC1'])
cor.test(pd[pd$Region=="ERC",'NeuN_pos'], pd[pd$Region=="ERC",'PC1'])
cor.test(pd[pd$Region=="HIPPO",'NeuN_pos'], pd[pd$Region=="HIPPO",'PC1'])

##
PC1_v_NeuN <- ggplot(data=pd, aes(x=PC1, y=NeuN_pos, col=Region) ) +
	 geom_point(size=5)+ 
	 labs(x=paste0("PC1 \n(",signif(varExplained[1],3) ,"% of variance)"),
		  y=paste0("Estimated Proportion NeuN+")) +   theme(legend.position = c(.12,.16), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2") + facet_wrap(~Region)  
		  
ggsave(PC1_v_NeuN, file='qc/SupplementalFigure_PC1_vs_Estimated_Prop_NeuNPlus_n377.pdf',height=10,width=12,useDingbats=FALSE )	 
############### Test correlation between estimated proportions via transcriptome and methylation ##########
pdMeth=pd
load('/dcl01/lieber/ajaffe/Steve/Alz/Paper/LIBD_AD_pd_hg38_withComp.Rdata',verbose=T)

pd$Neurons
pd$Dx <- plyr::revalue(pd$Dx, c('NC'='Control','AD'='Alzheimer') )
pd$Region <- plyr::revalue(pd$Region, c('CB'='CRB') )

summary(lm(Neurons~Dx+Region,data=pd))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="CRB",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="ERC",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="DLPFC",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="HIPPO",]))

pdf('plots/SupplementalFigure_RNAseq_estimated_proportion_of_neurons_by_dx_by_region.pdf',height=12,width=12,useDingbats=FALSE)
cell_type_dx_rnaseq = ggplot(data=pd, aes(x=Region,y=Neurons,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85),size=2) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Brain Region', y = "Estimated Proportion of Neurons") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
  theme(legend.position = 'none', 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )
print(cell_type_dx_rnaseq)
dev.off()

#################
pd_rna = pd 
pd_rna$Region = plyr::revalue(pd_rna$Region, c("CB"="CRB") )
pd_rna$Dx = plyr::revalue(pd_rna$Dx, c("NC"="Control","AD" = "Alzheimer") )
pd_rna$BrNum = paste0("Br", pd_rna$BRNum)

##
pdMeth$matchBrain = paste0(pdMeth$BrNum, "_", pdMeth$Region)
pd_rna$matchBrain = paste0(pd_rna$BrNum, "_", pd_rna$Region)

shared_samples = intersect(pdMeth$matchBrain,pd_rna$matchBrain )
pdMeth = pdMeth[pdMeth$matchBrain%in%shared_samples,]
pd_rna = pd_rna[pd_rna$matchBrain%in%shared_samples,]

### Dropping duplicated brains in RNAseq data
library(dplyr)
dropBrains = pd_rna[duplicated(pd_rna$matchBrain) | duplicated(pd_rna$matchBrain, fromLast=TRUE), ] 
dropBrains = dropBrains[order(dropBrains$BrNum, dropBrains$mitoRate),]
dropBrains = dropBrains[duplicated(dropBrains$matchBrain),'SAMPLE_ID']
pd_rna = pd_rna[!pd_rna$SAMPLE_ID %in% dropBrains, ]

### 
aa=match(shared_samples,pd_rna$matchBrain)
pd_rna=pd_rna[aa,]
bb=match(shared_samples,pdMeth$matchBrain)
pdMeth=pdMeth[bb,]
dat=cbind(pd_rna,pdMeth[,c('ES', 'NPC', 'DA_NEURON', 'NeuN_pos','NeuN_neg')])

cor.test(dat[,'Neurons'], dat[,'NeuN_pos'])
cor.test(dat[dat$Region=="CRB",'Neurons'], dat[dat$Region=="CRB",'NeuN_pos'])
cor.test(dat[dat$Region=="ERC",'Neurons'], dat[dat$Region=="ERC",'NeuN_pos'])
cor.test(dat[dat$Region=="HIPPO",'Neurons'], dat[dat$Region=="HIPPO",'NeuN_pos'])
cor.test(dat[dat$Region=="DLPFC",'Neurons'], dat[dat$Region=="DLPFC",'NeuN_pos'])

summary(lm(Neurons~NeuN_pos,data=dat[dat$Region=="CRB",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="ERC",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="DLPFC",]))
summary(lm(Neurons~Dx+totalAssignedGene+mitoRate+age+Sex+Race+RIN,data=pd[pd$Region=="HIPPO",]))

Neuron_v_NeuN <- ggplot(data=dat, aes(x=Neurons, y=NeuN_pos, col=Dx) ) +facet_wrap(~Region,nrow=4) +
	 geom_point(size=3)+ 
	 labs(x=paste0("RNA-seq Estimated Neuronal Proportion"),
		  y=paste0("DNAm Estimated Neuronal Proportion")) + scale_colour_brewer(palette="Dark2") +theme(legend.position='bottom', legend.title=element_blank()) + 
		  scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) 

ggsave(Neuron_v_NeuN, file='plots/SupplementalFigure_DNAm_vs_transcriptome_NeuN_Estimates_n264.pdf',height=20,width=18,useDingbats=FALSE )	 

Neuron_v_NeuN <- ggplot(data=dat, aes(x=Neurons, y=NeuN_pos, col=Dx) ) +
	 geom_point(size=3)+ 
	 labs(x=paste0("RNA-seq Estimated Neuronal Proportion"),
		  y=paste0("DNAm Estimated Neuronal Proportion")) +
scale_colour_manual(values=c("black","#ab1323" ) )
ggsave(Neuron_v_NeuN, file='plots/SupplementalFigure_DNAm_vs_transcriptome_NeuN_Estimates_n264_Dx.pdf',height=12,width=12,useDingbats=FALSE )	 


cor.test(dat[,'Neurons'],dat[,'NeuN_pos'])

corMat=cor( dat[,c('Neurons','OPC', 'Astrocytes','Oligodendrocytes','Microglia','Endothelial', 'ES','NPC','DA_NEURON','NeuN_pos','NeuN_neg')])
library(pheatmap)
pheatmap(corMat)


######### Checking aging information ####################
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper/')

## Aging analysis
load('rdas/aging_controlOnly_mergedStats.rda')

##
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')
load('rdas/merged_DMP_regionSpecific_caseControl_stats.rda')

tab = table(allRegion_mergedStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05,Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]<1e-3 & sign(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_logFC"]) == sign(allRegion_mergedStats[,'ALL_subset_mainEffect_logFC']))
fisher.test(tab)
tab

sig = allRegion_mergedStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05 & Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]<1e-3 & sign(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_logFC"]) == sign(allRegion_mergedStats[,'ALL_subset_mainEffect_logFC'])
allRegion_mergedStats[sig,'within10kb_geneSymbol_gencode_hg38']

aka = allRegion_mergedStats[sig,]
aka[aka$`within10kb_geneSymbol_gencode_hg38.1`%in% c("BIN1","ANK1","ANKRD30B","WDR81", "MYO1C","RHBDF2"),'Name']

Aging_controlOnly_mergedStats[aka[grep("BIN1",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value')]
Aging_controlOnly_mergedStats[aka[grep("ANK1",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']
Aging_controlOnly_mergedStats[aka[grep("WDR81",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']
Aging_controlOnly_mergedStats[aka[grep("MYO1C",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']
Aging_controlOnly_mergedStats[aka[grep("RHBDF2",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']

Aging_controlOnly_mergedStats[aka[grep("CSNK1G2",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']

Aging_controlOnly_mergedStats['cg19803550' , 'ALL_Control_Aging_P.Value']  #WDR81
Aging_controlOnly_mergedStats['cg05066959' , 'ALL_Control_Aging_P.Value']  #ANK1
Aging_controlOnly_mergedStats['cg14462670' , 'ALL_Control_Aging_P.Value']  #MYO1C

   #
  #

Aging_controlOnly_mergedStats[aka[grep("ANKRD30B",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']
Aging_controlOnly_mergedStats[aka[grep("DUSP22",aka$`within10kb_geneSymbol_gencode_hg38.1`),'Name'] , 'ALL_Control_Aging_P.Value']

