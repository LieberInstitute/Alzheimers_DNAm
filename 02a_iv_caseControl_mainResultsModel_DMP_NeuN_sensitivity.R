#qsub -l bluejay,mf=30G,h_vmem=40G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 02a_iv_caseControl_mainResultsModel_DMP_NeuN_sensitivity.R

#### Modeling methylation changes
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

###### All region analysis: Main models ######

###--- Subset ---###
## Main effect
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx + NeuN_pos + negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex,])

corfit <- duplicateCorrelation(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex])
fit <- lmFit(bVals[, subsetIndex], mod, block=pd$BrNum[subsetIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
subset_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2, genelist=goldsetSub)
table(subset_ALL_mainEffect$adj.P.Val<0.10)
colnames(subset_ALL_mainEffect)[25:30] <- paste0("ALL_subset_mainEffect_",colnames(subset_ALL_mainEffect)[25:30])

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
ord_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2)
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
col_order = c(colnames(allRegion_mergedStats)[1:24], 
			grep("_adj.P.Val",colnames(allRegion_mergedStats),value=T),	  
			grep("_logFC",colnames(allRegion_mergedStats),value=T),
			grep("_P.Value",colnames(allRegion_mergedStats),value=T),
			grep("_t$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_F$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(allRegion_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(allRegion_mergedStats),value=T,fixed=F) )
			
allRegion_mergedStats = allRegion_mergedStats[ ,col_order]
save(allRegion_mergedStats, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_mergedStats_DMP_analysis_dupCor_NeuN_sensitivity.rda' )		
	
####--- Full ---###
#
### Main effect
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd)
#fit <- lmFit(bVals, mod)
#fitEb <- eBayes(fit)
#full_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2)
#table(full_ALL_mainEffect$adj.P.Val<0.10)
#colnames(full_ALL_mainEffect) <- paste0("ALL_full_mainEffect_",colnames(full_ALL_mainEffect))
#
### Interaction effect
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region + Dx:Region, data = pd)
#fit <- lmFit(bVals, mod)
#fitEb <- eBayes(fit)
#
#full_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=11:13)
#table(full_ALL_interactionEffect$adj.P.Val<0.10)
#colnames(full_ALL_interactionEffect) <- paste0("ALL_full_interactionEffect_",colnames(full_ALL_interactionEffect))
#full_ALL_interactionEffect = full_ALL_interactionEffect[,-(1:3)]