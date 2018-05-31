library(limma)
library('bumphunter')
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(jaffelab)
source('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/00_brainseq_phase2_methylation_functions.R')				 
#
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')				  					  
					  
# drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
length(drop_hg38_unmappable)
bVals <- bVals[-drop_hg38_unmappable, ] 

# reorder goldset$Name
goldsetSub = goldset[match(rownames(bVals),goldset$Name),]
					  
##
#pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

### Get metastatistics for DMRs
load('rdas/merged_DMR.rda')

library(tidyr)
library(dplyr)
DMR_MetaStats = group_by(mergedDMR, Model, Region) %>% summarise(N_DMR_FWER_01 = sum (fwer <0.01), 
													  N_DMR_FWER_05 = sum (fwer <0.05 ), 
													  N_DMR_FWER_10 = sum (fwer <0.10 ),
													  DMR_Up_N =sum(value[fwer <0.10]>0), 
													  DMR_Down_N=sum(value[fwer <0.10]<0), 
													  DMR_Up_Percent = 100*DMR_Up_N/N_DMR_FWER_10, 
													  DMR_Median_Length = median(end[fwer <0.10]-start[fwer <0.10]),
													  DMR_Median_Area = median(area[fwer <0.10])													  
													  ) %>% as.data.frame()
DMR_MetaStats$DMR_Direction_P = apply(DMR_MetaStats[,c('DMR_Up_N', 'DMR_Down_N')] ,1, function(x){ try(chisq.test(x)$p.value) } )
DMR_MetaStats$DMR_Direction_P[grep("Error",DMR_MetaStats$DMR_Direction_P)] = NA
DMR_MetaStats$DMR_Direction_P <- as.numeric(DMR_MetaStats$DMR_Direction_P)
write.csv(DMR_MetaStats, 'csvs/DMR_summary_statistics.csv',row.names=F)

#### DMR Plots
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_annotation/GenomicState.Hsapiens.UCSC.hg38.knownGene.rda')

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
theTranscripts = annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene,codingOnly=TRUE)

an = annotateNearest(mergedDMR, theTranscripts)
mergedDMR$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
mergedDMR$nearestGeneDist = an$dist
mergedDMR$width=mergedDMR$end-mergedDMR$start
write.csv(mergedDMR[mergedDMR$fwer<0.05,],'csvs/SupplementalTable_Significant_DMRs.csv',row.names=F) 

##
subsetMain_DMR =  mergedDMR[mergedDMR$Model=="Subset Main", ]
sigDMR_subsetMain = mergedDMR[mergedDMR$cluster%in% subsetMain_DMR$cluster[subsetMain_DMR$fwer<0.05],]
sigDMR_subsetMain = sigDMR_subsetMain[order(sigDMR_subsetMain$cluster,sigDMR_subsetMain$Model,sigDMR_subsetMain$Region),]
write.csv(sigDMR_subsetMain,'csvs/SupplementalTable_sigDMRs_regionSpecificModels.csv',row.names=F) 

#######
cl = clusterMaker(as.character(goldsetSub$chr_hg38), goldsetSub$pos_hg38, maxGap=500)

dat= mergedDMR[mergedDMR$Model == "Subset Main", ] 
dat[order(dat$fwer),]

genes = matchGenes(dat[1:20,], theTranscripts)

subjectIndex = which(pd$DxOrdinal!="Alz Drop")
coi = paste0(pd$Region, ": ", pd$Dx)

mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region, data = pd[subjectIndex,])
cleaned_bVals= jaffelab::cleaningY( y = bVals[,subjectIndex, drop=FALSE], mod = mod, P = 2 )

### Various DMR plots
pdf('plots/SupplementalFigure_subset_allRegions_DUSP22.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[c(1),], p=cleaned_bVals, chr=goldsetSub$chr_hg38, pos=goldsetSub$pos_hg38, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()

pdf('plots/Figure_subset_allRegions_ANKRD30B.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[c(2),], p=cleaned_bVals, chr=goldsetSub$chr_hg38, pos=goldsetSub$pos_hg38, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()


pdf('plots/SupplementalFigure_subset_allRegions_JRK.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[c(3),], p=cleaned_bVals, chr=goldsetSub$chr_hg38, pos=goldsetSub$pos_hg38, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()

pdf('plots/SupplementalFigure_subset_allRegions_NAPRT.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[c(4),], p=cleaned_bVals, chr=goldsetSub$chr_hg38, pos=goldsetSub$pos_hg38, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()


pdf('plots/Figure_subset_allRegions_top_DMRs.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[1:20,], p=cleaned_bVals, chr=goldsetSub$chr_hg38, pos=goldsetSub$pos_hg38, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()

##
load('rdas/merged_DMR_NeuN_Sensitivity.rda',verbose=T)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
theTranscripts = annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene,codingOnly=TRUE)

an = annotateNearest(mergedDMR, theTranscripts)
mergedDMR$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
mergedDMR$nearestGeneDist = an$dist
mergedDMR$width=mergedDMR$end-mergedDMR$start
mergedDMR[mergedDMR$fwer<0.05,]