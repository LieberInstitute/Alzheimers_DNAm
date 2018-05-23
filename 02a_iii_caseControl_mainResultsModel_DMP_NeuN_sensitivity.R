#qsub -l bluejay,mf=30G,h_vmem=40G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 02a_iii_caseControl_mainResultsModel_DMP_NeuN_sensitivity.R

#### Modeling methylation changes
library(minfi)
library(limma)
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')
###

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
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx + NeuN_pos + negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex,])

corfit <- duplicateCorrelation(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
subset_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2, genelist=goldsetSub, confint=TRUE)
table(subset_ALL_mainEffect$adj.P.Val<0.10)
colnames(subset_ALL_mainEffect)[40:47] <- paste0("ALL_subset_mainEffect_",colnames(subset_ALL_mainEffect)[40:47])

## Interaction effect
mod <- model.matrix(~Dx + NeuN_pos + negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region + Dx:Region, data = pd[subsetIndex,])
fit <- lmFit(bVals[,subsetIndex], mod)
fitEb <- eBayes(fit)

subset_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=12:14)
table(subset_ALL_interactionEffect$adj.P.Val<0.10)
colnames(subset_ALL_interactionEffect) <- paste0("ALL_subset_interactionEffect_",colnames(subset_ALL_interactionEffect))
subset_ALL_interactionEffect = subset_ALL_interactionEffect[,-(1:3)]

###--- Ordinal ---###
## Main effect
mod <- model.matrix(~DxOrdinal + NeuN_pos + negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd)

corfit <- duplicateCorrelation(bVals, mod, block=pd$BrNum)
fit <- lmFit(bVals, mod, block=pd$BrNum, correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
ord_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ord_ALL_mainEffect$adj.P.Val<0.10)
colnames(ord_ALL_mainEffect) <- paste0("ALL_ordinal_mainEffect_",colnames(ord_ALL_mainEffect))

## Interaction effect
mod <- model.matrix(~DxOrdinal + NeuN_pos + negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region + DxOrdinal:Region, data = pd)

fit <- lmFit(bVals, mod)
fitEb <- eBayes(fit)

ord_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=12:14)
table(ord_ALL_interactionEffect$adj.P.Val<0.10)
colnames(ord_ALL_interactionEffect) <- paste0("ALL_ordinal_interactionEffect_",colnames(ord_ALL_interactionEffect))
ord_ALL_interactionEffect = ord_ALL_interactionEffect[,-(1:3)]

######### Merge
baseRownames = rownames(subset_ALL_mainEffect)

allRegion_mergedStats=
cbind(subset_ALL_mainEffect, 
	  subset_ALL_interactionEffect[baseRownames, ],
	  ord_ALL_mainEffect[baseRownames, ],
	  ord_ALL_interactionEffect[baseRownames, ] )
	  
### Reorder columns
col_order = c(colnames(allRegion_mergedStats)[1:39], 
			grep("_adj.P.Val",colnames(allRegion_mergedStats),value=T),	  
			grep("_logFC",colnames(allRegion_mergedStats),value=T),
			grep("_P.Value",colnames(allRegion_mergedStats),value=T),
			grep("_t$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_F$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(allRegion_mergedStats),value=T,fixed=F) )
			
allRegion_mergedStats = allRegion_mergedStats[ ,col_order]
save(allRegion_mergedStats, file='rdas/allRegion_mergedStats_DMP_analysis_dupCor_NeuN_sensitivity.rda' )		