library(limma)
library('bumphunter')
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(jaffelab)
source('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/00_brainseq_phase2_methylation_functions.R')				 
#
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')				  
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
					  
					  
# drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
length(drop_hg38_unmappable)
bVals <- bVals[-drop_hg38_unmappable, ] 

# reorder goldset$Name
goldsetSub = goldset[match(rownames(bVals),goldset$Name),]
					  
##
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )

##
pd$keepList[is.na(pd$keepList)] = TRUE
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
#pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

### Get metastatistics for DMRs
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMR.rda')

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
#write.csv(DMR_MetaStats, '/dcl01/lieber/ajaffe/Steve/Alz/csvs/DMR_summary_statistics.csv',row.names=F)

#### DMR Plots
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_annotation/GenomicState.Hsapiens.UCSC.hg38.knownGene.rda')

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
theTranscripts = annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene,codingOnly=TRUE)

an = annotateNearest(mergedDMR, theTranscripts)
mergedDMR$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
mergedDMR$nearestGeneDist = an$dist
mergedDMR$width=mergedDMR$end-mergedDMR$start
write.csv(mergedDMR[mergedDMR$fwer<0.10,],'/dcl01/lieber/ajaffe/Steve/Alz/csvs/Significant_DMRs.csv',row.names=F) 

##
subsetMain_DMR =  mergedDMR[mergedDMR$Model=="Subset Main", ]
sigDMR_subsetMain = mergedDMR[mergedDMR$cluster%in% subsetMain_DMR$cluster[subsetMain_DMR$fwer<0.1],]
sigDMR_subsetMain = sigDMR_subsetMain[order(sigDMR_subsetMain$cluster,sigDMR_subsetMain$Model,sigDMR_subsetMain$Region),]
#write.csv(sigDMR_subsetMain,'/dcl01/lieber/ajaffe/Steve/Alz/csvs/subsetMain_sigDMRs_acrossModels.csv',row.names=F) 

#######
cl = clusterMaker(as.character(goldsetSub$chr_hg38), goldsetSub$predictedPos, maxGap=500)

dat= mergedDMR[mergedDMR$Model == "Subset Main", ] 
dat[order(dat$fwer),]

genes = matchGenes(dat[1:20,], theTranscripts)

subjectIndex = which(pd$DxOrdinal!="Alz Drop")
coi = paste0(pd$Region, ": ", pd$Dx)

mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region, data = pd[subjectIndex,])
cleaned_bVals= jaffelab::cleaningY( y = bVals[,subjectIndex, drop=FALSE], mod = mod, P = 2 )


pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/subset_allRegions_top_DMRs.pdf',w=9,useDingbats=FALSE)
dmrPlot(regions=dat[1:20,], p=cleaned_bVals[], chr=goldsetSub$chr_hg38, pos=goldsetSub$predictedPos, cluster=cl, genomicState=GenomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, coi=pd$Dx[subjectIndex], genes=genes,  cols = c("black","#ab1323" ), Jitter=TRUE,linesSmooth=TRUE, build='hg38')
dev.off()
