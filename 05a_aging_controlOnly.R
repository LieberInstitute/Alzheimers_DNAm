#qsub -l bluejay,mf=50G,h_vmem=60G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 05a_aging_controlOnly.R

###### Modeling methylation changes
library(minfi)
library(limma)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
###
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  
goldsetSub = plyr::rename(goldsetSub, c('predictedPos'='pos_hg38','pos'='pos_hg19','chr'='chr_hg19') )
goldsetSub = goldsetSub[,c('chr_hg19','pos_hg19','chr_hg38','pos_hg38', intersect(colnames(goldsetSub), colnames(ann450kSub)) )]		
			  
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd$keepList[is.na(pd$keepList)] = TRUE
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))
			  
############# Controls only
get_agingStats_control = function(region) {
regionIndex=which(pd$Region==region & pd$Dx=="Control" )
mod <- model.matrix(~Age+ negControl_PC1 + negControl_PC2 +Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
stats <- topTable(fitEb, num=Inf, coef='Age', confint=TRUE,sort.by="none")
colnames(stats)<-paste0(region, "_Control_Aging_",colnames(stats))
return(stats)
}
controlAgingStats = lapply(unique(pd$Region), get_agingStats_control)		  
controlAgingStats = do.call("cbind",controlAgingStats)			  
colSums(controlAgingStats[,grep("adj.P.Val",colnames(controlAgingStats) )]<0.05)
			  
### main effect model
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Age + negControl_PC1 + negControl_PC2 +Sex + snpPC1 + Region, data = pd[subsetIndex,])

corfit <- duplicateCorrelation(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
ageMainControl <- topTable(fitEb, num=Inf, coef='Age',sort.by="none")
colnames(ageMainControl) <- paste0( "ALL_Control_Aging_", colnames(ageMainControl) )
table(ageMainControl$ALL_Control_Aging_adj.P.Val<0.05)

### interaction model
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Age + negControl_PC1 + negControl_PC2 +Sex + snpPC1 + Region + Age:Region, data = pd[subsetIndex,])

corfit <- duplicateCorrelation(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
ageIntControl <- topTable(fitEb, num=Inf, coef=10:12,sort.by="none")
colnames(ageIntControl) <- paste0( "INT_Control_Aging_", colnames(ageIntControl) )
table(ageIntControl$INT_Control_Aging_adj.P.Val<0.05 )

############## Merge all these tests ############
baseRownames = rownames(ageMainControl)
Aging_controlOnly_mergedStats=
cbind(ageMainControl[baseRownames, ],
	  ageIntControl[baseRownames, ],
	  controlAgingStats[baseRownames, ])
	  
### Reorder columns
col_order = c(grep("_adj.P.Val",colnames(Aging_controlOnly_mergedStats),value=T),	  
			grep("_logFC",colnames(Aging_controlOnly_mergedStats),value=T),
			grep("_P.Value",colnames(Aging_controlOnly_mergedStats),value=T),
			grep("_t$",colnames(Aging_controlOnly_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(Aging_controlOnly_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(Aging_controlOnly_mergedStats),value=T,fixed=F),
			grep("_CI",colnames(Aging_controlOnly_mergedStats),value=T,fixed=F) )
Aging_controlOnly_mergedStats=Aging_controlOnly_mergedStats[,col_order]			

save(Aging_controlOnly_mergedStats, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/aging_controlOnly_mergedStats.rda')	