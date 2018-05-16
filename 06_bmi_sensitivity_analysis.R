### BMI Sensitivity analysis

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
			  
## 
library(LIBDpheno)			  
pd = cbind(pd, bmi = toxicant[[1]][match(brnum(pd$BrNum),toxicant[[1]]$brnumerical),'bmi_calculated'])
			  
############# Controls only
get_bmiStats_control = function(region) {
regionIndex=which(pd$Region==region & pd$Dx=="Control" &!is.na(pd$bmi) )
mod <- model.matrix(~bmi+Age+ negControl_PC1 + negControl_PC2 +Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
stats <- topTable(fitEb, num=Inf, coef='bmi', confint=TRUE,sort.by="none")
colnames(stats)<-paste0(region, "_Control_BMI_",colnames(stats))
return(stats)
}
controlBMIStats = lapply(unique(pd$Region), get_bmiStats_control)		  
controlBMIStats = do.call("cbind",controlBMIStats)			  
colSums(controlBMIStats[,grep("adj.P.Val",colnames(controlBMIStats) )]<0.05)

toCheck=controlBMIStats[rowSums(controlBMIStats[,grep("adj.P.Val",colnames(controlBMIStats) )]<0.20)>0,]                                          