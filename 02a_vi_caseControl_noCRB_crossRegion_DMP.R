#qsub -l bluejay,mf=60G,h_vmem=60G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 02a_vi_caseControl_noCRB_crossRegion_DMP.R

#### Modeling methylation changes
library(minfi)
library(limma)
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  
			  
pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

###### All region analysis: Main models ######

###--- Subset ---###
## Main effect
subsetIndex=which(pd$keepList &pd$Region!="CRB" )
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + droplevels(Region), data = pd[subsetIndex,])

corfit <- duplicateCorrelation(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
subset_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2, genelist=goldsetSub, confint=TRUE)
table(subset_ALL_mainEffect$adj.P.Val<0.10)
colnames(subset_ALL_mainEffect)[40:47] <- paste0("Cortex_subset_mainEffect_",colnames(subset_ALL_mainEffect)[40:47])


###--- Ordinal ---###
## Main effect
subsetIndex=which(pd$Region!="CRB" )
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex])

corfit <- duplicateCorrelation(bVals[,subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[,subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
ord_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ord_ALL_mainEffect$adj.P.Val<0.10)
colnames(ord_ALL_mainEffect) <- paste0("Cortex_ordinal_mainEffect_",colnames(ord_ALL_mainEffect))


######### Merge
baseRownames = rownames(subset_ALL_mainEffect)

allRegion_mergedStats=
cbind(subset_ALL_mainEffect, 
	  ord_ALL_mainEffect[baseRownames, ] )
	  
### Reorder columns
col_order = c(colnames(allRegion_mergedStats)[1:39], 
			grep("_adj.P.Val",colnames(allRegion_mergedStats),value=T),	  
			grep("_logFC",colnames(allRegion_mergedStats),value=T),
			grep("_P.Value",colnames(allRegion_mergedStats),value=T),
			grep("_t$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_F$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(allRegion_mergedStats),value=T,fixed=F) )

cortexRegion_mergedStats = allRegion_mergedStats[ ,col_order]
save(cortexRegion_mergedStats, file='rdas/cortexRegions_mergedStats_DMP_analysis_dupCor.rda')