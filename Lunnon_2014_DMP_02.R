#qsub -l bluejay,mf=50G,h_vmem=60G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save Lunnon_2014_DMP_02.R

library(limma)

load('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/processed_bVals.rda')
bVals = processed_bVals

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

### drop hg38 unmappable
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  
goldsetSub = plyr::rename(goldsetSub, c('predictedPos'='pos_hg38','pos'='pos_hg19','chr'='chr_hg19') )
goldsetSub = goldsetSub[,c('chr_hg19','pos_hg19','chr_hg38','pos_hg38', intersect(colnames(goldsetSub), colnames(ann450kSub)) )]		

### drop sex probes
drop_sex_probes = which(rownames(bVals) %in% goldsetSub[as.character(goldsetSub$chr_hg38)%in%c("chrX","chrY"),'Name'])
length(drop_sex_probes)
bVals <- bVals[-drop_sex_probes, ] 					  

### drop cross reactive probes
xReactiveProbes <- read.csv(file="/users/ssemick/Tox_Expr/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
drop_xreactive <- which(rownames(bVals) %in% xReactiveProbes$TargetID)
length(drop_xreactive)
bVals <- bVals[-drop_xreactive,] 

###### Regional analysis disease ######

####---DLPFC---####
regionIndex=which(pd$Region=="PFC" & !is.na(pd$Dx) )
mod <- model.matrix(~Dx, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
DxStats_PFC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(DxStats_PFC$adj.P.Val<0.10)
colnames(DxStats_PFC) <- paste0("PFC_Dx_", colnames(DxStats_PFC) )

####---STG---####
regionIndex=which(pd$Region=="STG" & !is.na(pd$Dx))
mod <- model.matrix(~Dx, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
DxStats_STG <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(DxStats_STG$adj.P.Val<0.10)
colnames(DxStats_STG) <- paste0("STG_Dx_",colnames(DxStats_STG))

####---ERC---####
regionIndex=which(pd$Region=="ERC" & !is.na(pd$Dx))
mod <- model.matrix(~Dx, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
DxStats_ERC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(DxStats_ERC$adj.P.Val<0.10)
colnames(DxStats_ERC) <- paste0("ERC_Dx_",colnames(DxStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB" & !is.na(pd$Dx))
mod <- model.matrix(~Dx, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
DxStats_CRB <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(DxStats_CRB$adj.P.Val<0.10)
colnames(DxStats_CRB) <- paste0("CRB_Dx_",colnames(DxStats_CRB))

####---ALL Regions ---####

## Main effect
regionIndex=which(!is.na(pd$Dx)& pd$Region!="blood")
mod <- model.matrix(~Dx +Region, data = pd[regionIndex,])

corfit <- duplicateCorrelation(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex])
fit <- lmFit(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
DxStats_ALL_main <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(DxStats_ALL_main$adj.P.Val<0.05)
colnames(DxStats_ALL_main) <- paste0("ALL_main_Dx_",colnames(DxStats_ALL_main))

## Interaction effect
regionIndex=which(!is.na(pd$Dx) & pd$Region!="blood")
mod <- model.matrix(~Dx +Region + Region:Dx, data = pd[regionIndex,])

corfit <- duplicateCorrelation(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex])
fit <- lmFit(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
DxStats_ALL_int <- topTable(fitEb, num=Inf, coef=6:8, confint=TRUE)
table(DxStats_ALL_int$adj.P.Val<0.10)
colnames(DxStats_ALL_int) <- paste0("ALL_int_Dx_",colnames(DxStats_ALL_int))

###### Regional analysis braak_stage ######

####---PFC---####
regionIndex=which(pd$Region=="PFC" & !is.na(pd$braak_stage) )
mod <- model.matrix(~braak_stage, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
BraakStats_PFC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(BraakStats_PFC$adj.P.Val<0.10)
colnames(BraakStats_PFC) <- paste0("PFC_Braak_", colnames(BraakStats_PFC) )

####---STG---####
regionIndex=which(pd$Region=="STG" & !is.na(pd$braak_stage))
mod <- model.matrix(~braak_stage, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
BraakStats_STG <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(BraakStats_STG$adj.P.Val<0.10)
colnames(BraakStats_STG) <- paste0("STG_Braak_",colnames(BraakStats_STG))

####---ERC---####
regionIndex=which(pd$Region=="ERC" & !is.na(pd$braak_stage))
mod <- model.matrix(~braak_stage, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
BraakStats_ERC <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(BraakStats_ERC$adj.P.Val<0.10)
colnames(BraakStats_ERC) <- paste0("ERC_Braak_",colnames(BraakStats_ERC))

####---CRB---####
regionIndex=which(pd$Region=="CRB" & !is.na(pd$braak_stage))
mod <- model.matrix(~braak_stage, data = pd[regionIndex,])
fit <- lmFit(bVals[,regionIndex], mod)
fitEb <- eBayes(fit)
BraakStats_CRB <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(BraakStats_CRB$adj.P.Val<0.10)
colnames(BraakStats_CRB) <- paste0("CRB_Braak_",colnames(BraakStats_CRB))

####---ALL Regions ---####

## Main effect
regionIndex=which(!is.na(pd$braak_stage)& pd$Region!="blood")
mod <- model.matrix(~braak_stage +Region, data = pd[regionIndex,])

corfit <- duplicateCorrelation(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex])
fit <- lmFit(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
BraakStats_ALL_main <- topTable(fitEb, num=Inf, coef=2, confint=TRUE)
table(BraakStats_ALL_main$adj.P.Val<0.10)
colnames(BraakStats_ALL_main) <- paste0("ALL_main_Braak_",colnames(BraakStats_ALL_main))

## Interaction effect
regionIndex=which(!is.na(pd$braak_stage) & pd$Region!="blood")
mod <- model.matrix(~braak_stage +Region + Region:braak_stage, data = pd[regionIndex,])

corfit <- duplicateCorrelation(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex])
fit <- lmFit(bVals[, regionIndex], mod, block=pd$SampleID[regionIndex], correlation = corfit$consensus.correlation)

fitEb <- eBayes(fit)
BraakStats_ALL_int <- topTable(fitEb, num=Inf, coef=6:8, confint=TRUE)
table(BraakStats_ALL_int$adj.P.Val<0.10)
colnames(BraakStats_ALL_int) <- paste0("ALL_int_Braak_",colnames(BraakStats_ALL_int))

############## Merge all these tests ############
baseRownames = rownames(DxStats_PFC)
Lunnon_mergedStats=
cbind(DxStats_PFC, 
	  DxStats_STG[baseRownames, ],
	  DxStats_ERC[baseRownames, ],
	  DxStats_CRB[baseRownames, ],
	  DxStats_ALL_main[baseRownames, ],
	  DxStats_ALL_int[baseRownames, ],
	  BraakStats_PFC[baseRownames, ],
	  BraakStats_STG[baseRownames, ],
	  BraakStats_ERC[baseRownames, ],
	  BraakStats_CRB[baseRownames, ],
	  BraakStats_ALL_main[baseRownames, ],
	  BraakStats_ALL_int[baseRownames, ])	  
### Reorder columns
col_order = c( 
			grep("_adj.P.Val",colnames(Lunnon_mergedStats),value=T),	  
			grep("_logFC",colnames(Lunnon_mergedStats),value=T),
			grep("_P.Value",colnames(Lunnon_mergedStats),value=T),
			grep("_t$",colnames(Lunnon_mergedStats),value=T,fixed=F),
			grep("_F$",colnames(Lunnon_mergedStats),value=T,fixed=F),
			grep("_AveExpr$",colnames(Lunnon_mergedStats),value=T,fixed=F),
			grep("_B$",colnames(Lunnon_mergedStats),value=T,fixed=F),
			grep("_CI",colnames(Lunnon_mergedStats),value=T,fixed=F)			)
Lunnon_mergedStats=Lunnon_mergedStats[,col_order]			

#### define standard error as (R.CI-L.CI)/3.92 <- 1.96*2
ci = Lunnon_mergedStats[,grep("CI",colnames(Lunnon_mergedStats),v=T)]
split_col<-seq(1,ncol(ci),by=2)
se = mapply(function(x,y) { (y-x) / (qnorm(.975, mean = 0, sd = 1, log = FALSE) *2)  } ,ci[,split_col],ci[,-split_col] )
colnames(se) <- gsub("CI.L","SE",colnames( se) )
Lunnon_mergedStats= cbind(Lunnon_mergedStats, se)
##########

save(Lunnon_mergedStats, file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_Lunnon_mergedStats.rda')	
