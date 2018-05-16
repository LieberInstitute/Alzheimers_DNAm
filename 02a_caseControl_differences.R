#### Modeling methylation changes
library(minfi)
library(limma)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
#
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
					  
###### Full within region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
fullStats_DLPFC <- topTable(fitEb, num=Inf, coef=2, genelist=ann450kSub)
table(fullStats_DLPFC$adj.P.Val<0.05)
colnames(fullStats_DLPFC)[23:28] <- paste0("DLPFC_full_",colnames(fullStats_DLPFC)[23:28])

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
fullStats_HIPPO <- topTable(fitEb, num=Inf, coef=2)
table(fullStats_HIPPO$adj.P.Val<0.05)
colnames(fullStats_HIPPO) <- paste0("HIPPO_full_",colnames(fullStats_HIPPO))

####---ERC---####
regionIndex=which(pd$Region=="ERC")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
fullStats_ERC <- topTable(fitEb, num=Inf, coef=2)
table(fullStats_ERC$adj.P.Val<0.05)
colnames(fullStats_ERC) <- paste0("ERC_full_",colnames(fullStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
fullStats_CRB <- topTable(fitEb, num=Inf, coef=2)
table(fullStats_CRB$adj.P.Val<0.05)
colnames(fullStats_CRB) <- paste0("CRB_full_",colnames(fullStats_CRB))

###### Keep within region analysis ######
pd$keepList[is.na(pd$keepList)] = TRUE
####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_DLPFC <- topTable(fitEb, num=Inf, coef=2)
table(subStats_DLPFC$adj.P.Val<0.05)
colnames(subStats_DLPFC) <- paste0("DLPFC_subset_",colnames(subStats_DLPFC))

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO" & pd$keepList )
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_HIPPO <- topTable(fitEb, num=Inf, coef=2)
table(subStats_HIPPO$adj.P.Val<0.05)
colnames(subStats_HIPPO) <- paste0("HIPPO_subset_",colnames(subStats_HIPPO))

####---ERC---####
regionIndex=which(pd$Region=="ERC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_ERC <- topTable(fitEb, num=Inf, coef=2)
table(subStats_ERC$adj.P.Val<0.05)
colnames(subStats_ERC) <- paste0("ERC_subset_",colnames(subStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
subStats_CRB <- topTable(fitEb, num=Inf, coef=2)
table(subStats_CRB$adj.P.Val<0.05)
colnames(subStats_CRB) <- paste0("CRB_subset_",colnames(subStats_CRB))

###### Ordinal model region analysis ######

pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_DLPFC <- topTable(fitEb, num=Inf, coef=2)
table(ordStats_DLPFC$adj.P.Val<0.05)
colnames(ordStats_DLPFC) <- paste0("DLPFC_ordinal_",colnames(ordStats_DLPFC))

####---Hippo---####
regionIndex=which(pd$Region=="HIPPO")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_HIPPO <- topTable(fitEb, num=Inf, coef=2)
table(ordStats_HIPPO$adj.P.Val<0.05)
colnames(ordStats_HIPPO) <- paste0("HIPPO_ordinal_",colnames(ordStats_HIPPO))

####---ERC---####
regionIndex=which(pd$Region=="ERC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_ERC <- topTable(fitEb, num=Inf, coef=2)
table(ordStats_ERC$adj.P.Val<0.05)
colnames(ordStats_ERC) <- paste0("ERC_ordinal_",colnames(ordStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
ordStats_CRB <- topTable(fitEb, num=Inf, coef=2)
table(ordStats_CRB$adj.P.Val<0.05)
colnames(ordStats_CRB) <- paste0("CRB_ordinal_",colnames(ordStats_CRB))

###### All region analysis ######

###--- Full ---###

## Main effect
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region, data = pd)
fit <- lmFit(bVals, mod)
fitEb <- eBayes(fit)
full_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2)
table(full_ALL_mainEffect$adj.P.Val<0.05)
colnames(full_ALL_mainEffect) <- paste0("ALL_full_mainEffect_",colnames(full_ALL_mainEffect))

## Interaction effect
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region + Dx:Region, data = pd)
fit <- lmFit(bVals, mod)
fitEb <- eBayes(fit)

full_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=11:13)
table(full_ALL_interactionEffect$adj.P.Val<0.05)
colnames(full_ALL_interactionEffect) <- paste0("ALL_full_interactionEffect_",colnames(full_ALL_interactionEffect))
full_ALL_interactionEffect = full_ALL_interactionEffect[,-(1:3)]

###--- Subset ---###
## Main effect
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region, data = pd[subsetIndex,])
fit <- lmFit(bVals[,subsetIndex], mod)
fitEb <- eBayes(fit)

subset_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2)
table(subset_ALL_mainEffect$adj.P.Val<0.05)
colnames(subset_ALL_mainEffect) <- paste0("ALL_subset_mainEffect_",colnames(subset_ALL_mainEffect))

## Interaction effect
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region + Dx:Region, data = pd[subsetIndex,])
fit <- lmFit(bVals[,subsetIndex], mod)
fitEb <- eBayes(fit)

subset_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=11:13)
table(subset_ALL_interactionEffect$adj.P.Val<0.05)
colnames(subset_ALL_interactionEffect) <- paste0("ALL_subset_interactionEffect_",colnames(subset_ALL_interactionEffect))
subset_ALL_interactionEffect = subset_ALL_interactionEffect[,-(1:3)]

###--- Ordinal ---###
## Main effect
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region, data = pd)
fit <- lmFit(bVals, mod)
fitEb <- eBayes(fit)
ord_ALL_mainEffect <- topTable(fitEb, num=Inf, coef=2)
table(ord_ALL_mainEffect$adj.P.Val<0.05)
colnames(ord_ALL_mainEffect) <- paste0("ALL_ordinal_mainEffect_",colnames(ord_ALL_mainEffect))

## Interaction effect
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + Race + Region + DxOrdinal:Region, data = pd)
fit <- lmFit(bVals, mod)
fitEb <- eBayes(fit)

ord_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=11:13)
table(ord_ALL_interactionEffect$adj.P.Val<0.05)
colnames(ord_ALL_interactionEffect) <- paste0("ALL_ordinal_interactionEffect_",colnames(ord_ALL_interactionEffect))
ord_ALL_interactionEffect = ord_ALL_interactionEffect[,-(1:3)]

############## Merge all these tests ############
baseRownames = rownames(fullStats_DLPFC)
mergedStats=
cbind(fullStats_DLPFC, 
	  fullStats_HIPPO[baseRownames, ],
	  fullStats_ERC[baseRownames, ],
	  fullStats_CRB[baseRownames, ],
	  subStats_DLPFC[baseRownames, ],
	  subStats_HIPPO[baseRownames, ],
	  subStats_ERC[baseRownames, ],
	  subStats_CRB[baseRownames, ],
	  ordStats_DLPFC[baseRownames, ],
	  ordStats_HIPPO[baseRownames, ],
	  ordStats_ERC[baseRownames, ],
	  ordStats_CRB[baseRownames, ],	  
	  full_ALL_mainEffect[baseRownames, ],
	  full_ALL_interactionEffect[baseRownames, ],
	  subset_ALL_mainEffect[baseRownames, ],
	  subset_ALL_interactionEffect[baseRownames, ],
	  ord_ALL_mainEffect[baseRownames, ],
	  ord_ALL_interactionEffect[baseRownames, ])
	  
### Reorder columns
col_order = c(colnames(mergedStats)[1:22], 
			grep("_adj.P.Val",colnames(mergedStats),value=T),	  
			grep("_logFC",colnames(mergedStats),value=T),
			grep("_P.Value",colnames(mergedStats),value=T),
			grep("_t$",colnames(mergedStats),value=T,fixed=F),
			grep("_F$",colnames(mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(mergedStats),value=T,fixed=F),
			grep("_B$",colnames(mergedStats),value=T,fixed=F) )
mergedStats=mergedStats[,col_order]			


save(mergedStats, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMC_case_control_stats.rda')	

tests = grep("adj.P.Val", colnames(mergedStats),value=TRUE )
table(rowSums(mergedStats[,tests] <0.10)>0)
write.csv(mergedStats[rowSums(mergedStats[,tests] <0.05)>0,],file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/merged_DMC_case_control_stats_FDR10.csv',row.names=FALSE)	
write.csv(mergedStats[rowSums(mergedStats[,tests] <0.05)>0,],file=gzfile('/dcl01/lieber/ajaffe/Steve/Alz/csvs/merged_DMC_case_control_stats_FDR10.csv.gz'),row.names=FALSE)	