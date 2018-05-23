#### Modeling methylation changes
library(minfi)
library(limma)

setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
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


###### Keep within region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_DLPFC <- topTable(fitEb, num=Inf, coef=2,genelist=goldsetSub, confint=TRUE)
table(subStats_DLPFC$adj.P.Val<0.10)
colnames(subStats_DLPFC)[40:47] <- paste0("DLPFC_subset_NoAdj_",colnames(subStats_DLPFC)[40:47])

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO" & pd$keepList )
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_HIPPO <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_HIPPO$adj.P.Val<0.10)
colnames(subStats_HIPPO) <- paste0("HIPPO_subset_NoAdj_",colnames(subStats_HIPPO))

####---ERC---####
regionIndex=which(pd$Region=="ERC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_ERC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_ERC$adj.P.Val<0.10)
colnames(subStats_ERC) <- paste0("ERC_subset_NoAdj_",colnames(subStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_CRB <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_CRB$adj.P.Val<0.10)
colnames(subStats_CRB) <- paste0("CRB_subset_NoAdj_",colnames(subStats_CRB))

###### Ordinal model region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_DLPFC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_DLPFC$adj.P.Val<0.10)
colnames(ordStats_DLPFC) <- paste0("DLPFC_ordinal_NoAdj_",colnames(ordStats_DLPFC))

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_HIPPO <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_HIPPO$adj.P.Val<0.10)
colnames(ordStats_HIPPO) <- paste0("HIPPO_ordinal_NoAdj_",colnames(ordStats_HIPPO))

####---ERC---####
regionIndex=which(pd$Region=="ERC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_ERC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_ERC$adj.P.Val<0.10)
colnames(ordStats_ERC) <- paste0("ERC_ordinal_NoAdj_",colnames(ordStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_CRB <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_CRB$adj.P.Val<0.10)
colnames(ordStats_CRB) <- paste0("CRB_ordinal_NoAdj_",colnames(ordStats_CRB))

####---Subset Sensitivity---####

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1+NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_DLPFC_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_DLPFC_adj$adj.P.Val<0.10)
colnames(subStats_DLPFC_adj) <- paste0("DLPFC_subset_cellTypeAdj_",colnames(subStats_DLPFC_adj))

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO" & pd$keepList )
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1+NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_HIPPO_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_HIPPO_adj$adj.P.Val<0.10)
colnames(subStats_HIPPO_adj) <- paste0("HIPPO_subset_cellTypeAdj_",colnames(subStats_HIPPO_adj))

####---ERC---####
regionIndex=which(pd$Region=="ERC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_ERC_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_ERC_adj$adj.P.Val<0.10)
colnames(subStats_ERC_adj) <- paste0("ERC_subset_cellTypeAdj_",colnames(subStats_ERC_adj))

####---CRB---####
regionIndex=which(pd$Region=="CRB" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_CRB_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(subStats_CRB_adj$adj.P.Val<0.10)
colnames(subStats_CRB_adj) <- paste0("CRB_subset_cellTypeAdj_",colnames(subStats_CRB_adj))

###### Ordinal model region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_DLPFC_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_DLPFC_adj$adj.P.Val<0.10)
colnames(ordStats_DLPFC_adj) <- paste0("DLPFC_ordinal_cellTypeAdj_",colnames(ordStats_DLPFC_adj))

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_HIPPO_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_HIPPO_adj$adj.P.Val<0.10)
colnames(ordStats_HIPPO_adj) <- paste0("HIPPO_ordinal_cellTypeAdj_",colnames(ordStats_HIPPO_adj))

####---ERC---####
regionIndex=which(pd$Region=="ERC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_ERC_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_ERC_adj$adj.P.Val<0.10)
colnames(ordStats_ERC_adj) <- paste0("ERC_ordinal_cellTypeAdj_",colnames(ordStats_ERC_adj))

####---CRB---####
regionIndex=which(pd$Region=="CRB")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + NeuN_pos, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_CRB_adj <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(ordStats_CRB_adj$adj.P.Val<0.10)
colnames(ordStats_CRB_adj) <- paste0("CRB_ordinal_cellTypeAdj_",colnames(ordStats_CRB_adj))

############## Merge all these tests ############
baseRownames = rownames(subStats_DLPFC)
regionSpecific_mergedStats=
cbind(subStats_DLPFC[baseRownames, ],
	  subStats_HIPPO[baseRownames, ],
	  subStats_ERC[baseRownames, ],
	  subStats_CRB[baseRownames, ],
	  subStats_DLPFC_adj[baseRownames, ],
	  subStats_HIPPO_adj[baseRownames, ],
	  subStats_ERC_adj[baseRownames, ],
	  subStats_CRB_adj[baseRownames, ],
	  ordStats_DLPFC[baseRownames, ],
	  ordStats_HIPPO[baseRownames, ],
	  ordStats_ERC[baseRownames, ],
	  ordStats_CRB[baseRownames, ],
	  ordStats_DLPFC_adj[baseRownames, ],
	  ordStats_HIPPO_adj[baseRownames, ],
	  ordStats_ERC_adj[baseRownames, ],
	  ordStats_CRB_adj[baseRownames, ])
	  
### Reorder columns
col_order = c(colnames(regionSpecific_mergedStats)[1:39], 
			grep("_adj.P.Val",colnames(regionSpecific_mergedStats),value=T),	  
			grep("_logFC",colnames(regionSpecific_mergedStats),value=T),
			grep("_P.Value",colnames(regionSpecific_mergedStats),value=T),
			grep("_t$",colnames(regionSpecific_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(regionSpecific_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(regionSpecific_mergedStats),value=T,fixed=F),
			grep("_CI",colnames(regionSpecific_mergedStats),value=T,fixed=F) )
regionSpecific_mergedStats=regionSpecific_mergedStats[,col_order]			

#### define standard error as (R.CI-L.CI)/3.92 <- 1.96*2
ci = regionSpecific_mergedStats[,grep("CI",colnames(regionSpecific_mergedStats),v=T)]
split_col<-seq(1,ncol(ci),by=2)
se = mapply(function(x,y) { (y-x) / (qnorm(.975, mean = 0, sd = 1, log = FALSE) *2)  } ,ci[,split_col],ci[,-split_col] )
colnames(se) <- gsub("CI.L","SE",colnames( se) )
regionSpecific_mergedStats= cbind(regionSpecific_mergedStats, se)

save(regionSpecific_mergedStats, file='rdas/merged_DMP_regionSpecific_caseControl_stats.rda')	