#qsub -l mf=40G,h_vmem=40G,h_fsize=200G,h_stack=256M -cwd -pe local 4 -R y -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 03a_caseControl_DMR_analysis.R

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

##
pd$keepList[is.na(pd$keepList)] = TRUE
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

### subset probes ###
# drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
length(drop_hg38_unmappable)
bVals <- bVals[-drop_hg38_unmappable, ] 

# reorder goldset$Name
goldsetSub = goldset[match(rownames(bVals),goldset$Name),]


#### DMR Analysis
library('doParallel')
registerDoParallel(cores = 4)
library('bumphunter')
#library('TxDb.Hsapiens.UCSC.hg19.knownGene')
#genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

###### Subset within region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

subsetDMRs_DLPFC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_DLPFC = subsetDMRs_DLPFC$table
subsetDMRs_DLPFC$Model = "Subset"
subsetDMRs_DLPFC$Region = "DLPFC"

####---HIPPO---####
regionIndex=which(pd$Region=="HIPPO"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

subsetDMRs_HIPPO <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_HIPPO = subsetDMRs_HIPPO$table
subsetDMRs_HIPPO$Model = "Subset"
subsetDMRs_HIPPO$Region = "HIPPO"

####---ERC---####
regionIndex=which(pd$Region=="ERC"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

subsetDMRs_ERC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_ERC = subsetDMRs_ERC$table
subsetDMRs_ERC$Model = "Subset"
subsetDMRs_ERC$Region = "ERC"


####---Cerebellum---####
regionIndex=which(pd$Region=="CRB"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

subsetDMRs_CRB <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_CRB = subsetDMRs_CRB$table
subsetDMRs_CRB$Model = "Subset"
subsetDMRs_CRB$Region = "CRB"
#####
regionalSubset = rbind(subsetDMRs_DLPFC, subsetDMRs_HIPPO, subsetDMRs_ERC, subsetDMRs_CRB)
save(regionalSubset, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/regionalSubset_DMR.rda')
###### Ordinal within region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

ordDMRs_DLPFC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
ordDMRs_DLPFC = ordDMRs_DLPFC$table
ordDMRs_DLPFC$Model = "Ordinal"
ordDMRs_DLPFC$Region = "DLPFC"

####---HIPPO---####
regionIndex=which(pd$Region=="HIPPO")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

ordDMRs_HIPPO <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
ordDMRs_HIPPO = ordDMRs_HIPPO$table
ordDMRs_HIPPO$Model = "Ordinal"
ordDMRs_HIPPO$Region = "HIPPO"

####---ERC---####
regionIndex=which(pd$Region=="ERC")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

ordDMRs_ERC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
ordDMRs_ERC = ordDMRs_ERC$table
ordDMRs_ERC$Model = "Ordinal"
ordDMRs_ERC$Region = "ERC"


####---Cerebellum---####
regionIndex=which(pd$Region=="CRB")
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])

ordDMRs_CRB <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
ordDMRs_CRB = ordDMRs_CRB$table
ordDMRs_CRB$Model = "Ordinal"
ordDMRs_CRB$Region = "CRB"

regionalOrdinal = rbind(ordDMRs_DLPFC, ordDMRs_HIPPO, ordDMRs_ERC, ordDMRs_CRB)
save(regionalOrdinal, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/regionalOrdinal_DMR.rda')
###### All region analysis ######

###--- Subset ---###

## Main effect
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex,])
subsetMainDMRs_ALL <- bumphunter( bVals[,subsetIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetMainDMRs_ALL = subsetMainDMRs_ALL$table
subsetMainDMRs_ALL$Model = "Subset Main"
subsetMainDMRs_ALL$Region = "ALL"
###--- Ordinal ---###

## Main effect
mod <- model.matrix(~DxOrdinal+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd)
ordinalMainDMRs_ALL <- bumphunter( bVals, 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$predictedPos,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
ordinalMainDMRs_ALL = ordinalMainDMRs_ALL$table
ordinalMainDMRs_ALL$Model = "Ordinal Main"
ordinalMainDMRs_ALL$Region = "ALL"

#####
allRegions = rbind(subsetMainDMRs_ALL, ordinalMainDMRs_ALL)
save(allRegions, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegions_DMR.rda')

#####
mergedDMR = rbind(regionalSubset,regionalOrdinal, allRegions )
save(mergedDMR, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMR.rda')

#########===== Full Analyses: No longer used =====#########

###--- Full ---###

### Main effect
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd)
#fullMainDMRs_ALL <- bumphunter( bVals, 
#							  design = mod, 
#							  cutoff = 0.05, 
#							  B=1000, 
#							  coef = 2,
#							  chr = goldsetSub$chr_hg38,
#							  pos = goldsetSub$predictedPos,
#							  nullMethod = 'bootstrap',
#							  smooth=TRUE)
#fullMainDMRs_ALL = fullMainDMRs_ALL$table
#fullMainDMRs_ALL$Model = "Full Main"
#fullMainDMRs_ALL$Region = "ALL"


###### Full within region analysis ######

#####---DLPFC---####
#regionIndex=which(pd$Region=="DLPFC")
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
#
#fullDMRs_DLPFC <- bumphunter( bVals[,regionIndex], 
#							  design = mod, 
#							  cutoff = 0.05, 
#							  B=1000, 
#							  coef = 2,
#							  chr = goldsetSub$chr_hg38,
#							  pos = goldsetSub$predictedPos,
#							  nullMethod = 'bootstrap',
#							  smooth=TRUE)
#fullDMRs_DLPFC = fullDMRs_DLPFC$table
#fullDMRs_DLPFC$Model = "Full"
#fullDMRs_DLPFC$Region = "DLPFC"
#
#####---HIPPO---####
#regionIndex=which(pd$Region=="HIPPO")
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
#
#fullDMRs_HIPPO <- bumphunter( bVals[,regionIndex], 
#							  design = mod, 
#							  cutoff = 0.05, 
#							  B=1000, 
#							  coef = 2,
#							  chr = goldsetSub$chr_hg38,
#							  pos = goldsetSub$predictedPos,
#							  nullMethod = 'bootstrap',
#							  smooth=TRUE)
#fullDMRs_HIPPO = fullDMRs_HIPPO$table
#fullDMRs_HIPPO$Model = "Full"
#fullDMRs_HIPPO$Region = "HIPPO"
#
#####---ERC---####
#regionIndex=which(pd$Region=="ERC")
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
#
#fullDMRs_ERC <- bumphunter( bVals[,regionIndex], 
#							  design = mod, 
#							  cutoff = 0.05, 
#							  B=1000, 
#							  coef = 2,
#							  chr = goldsetSub$chr_hg38,
#							  pos = goldsetSub$predictedPos,
#							  nullMethod = 'bootstrap',
#							  smooth=TRUE)
#fullDMRs_ERC = fullDMRs_ERC$table
#fullDMRs_ERC$Model = "Full"
#fullDMRs_ERC$Region = "ERC"
#
#
#####---Cerebellum---####
#regionIndex=which(pd$Region=="CRB")
#mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[regionIndex,])
#
#fullDMRs_CRB <- bumphunter( bVals[,regionIndex], 
#							  design = mod, 
#							  cutoff = 0.05, 
#							  B=1000, 
#							  coef = 2,
#							  chr = goldsetSub$chr_hg38,
#							  pos = goldsetSub$predictedPos,
#							  nullMethod = 'bootstrap',
#							  smooth=TRUE)
#fullDMRs_CRB = fullDMRs_CRB$table
#fullDMRs_CRB$Model = "Full"
#fullDMRs_CRB$Region = "CRB"
#
######
#regionalFull = rbind(fullDMRs_DLPFC,fullDMRs_HIPPO,fullDMRs_ERC,fullDMRs_CRB)
#save(regionalFull, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/regionalFull_DMR.rda')