#qsub -l mf=50G,h_vmem=50G,h_fsize=200G,h_stack=256M -cwd -pe local 4 -R y -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 03c_caseControl_DMR_analysis_NeuN_sensitivity.R

#### Modeling methylation changes
library(minfi)
library(limma)
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')
#

pd$DxOrdinal= as.numeric(factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") ))

### subset probes ###
# drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
length(drop_hg38_unmappable)
bVals <- bVals[-drop_hg38_unmappable, ] 

# reorder goldset$Name
goldsetSub = goldset[match(rownames(bVals),goldset$Name),]

#### DMR Analysis
library('doParallel')
registerDoParallel(cores = 4)
library('bumphunter')

###### Subset within region analysis ######

####---DLPFC---####
regionIndex=which(pd$Region=="DLPFC" & pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 +NeuN_pos, data = pd[regionIndex,])

subsetDMRs_DLPFC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$pos_hg38,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_DLPFC = subsetDMRs_DLPFC$table
subsetDMRs_DLPFC$Model = "Subset"
subsetDMRs_DLPFC$Region = "DLPFC"

####---HIPPO---####
regionIndex=which(pd$Region=="HIPPO"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1+NeuN_pos, data = pd[regionIndex,])

subsetDMRs_HIPPO <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$pos_hg38,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_HIPPO = subsetDMRs_HIPPO$table
subsetDMRs_HIPPO$Model = "Subset"
subsetDMRs_HIPPO$Region = "HIPPO"

####---ERC---####
regionIndex=which(pd$Region=="ERC"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1+NeuN_pos, data = pd[regionIndex,])

subsetDMRs_ERC <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$pos_hg38,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_ERC = subsetDMRs_ERC$table
subsetDMRs_ERC$Model = "Subset"
subsetDMRs_ERC$Region = "ERC"


####---Cerebellum---####
regionIndex=which(pd$Region=="CRB"& pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1+NeuN_pos, data = pd[regionIndex,])

subsetDMRs_CRB <- bumphunter( bVals[,regionIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$pos_hg38,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetDMRs_CRB = subsetDMRs_CRB$table
subsetDMRs_CRB$Model = "Subset"
subsetDMRs_CRB$Region = "CRB"
#####
regionalSubset = rbind(subsetDMRs_DLPFC, subsetDMRs_HIPPO, subsetDMRs_ERC, subsetDMRs_CRB)
save(regionalSubset, file='rdas/regionalSubset_DMR_NeuN_Sensitivity.rda')

###### All region analysis ######

###--- Subset ---###

## Main effect
subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region+NeuN_pos, data = pd[subsetIndex,])
subsetMainDMRs_ALL <- bumphunter( bVals[,subsetIndex], 
							  design = mod, 
							  cutoff = 0.05, 
							  B=1000, 
							  coef = 2,
							  chr = goldsetSub$chr_hg38,
							  pos = goldsetSub$pos_hg38,
							  nullMethod = 'bootstrap',
							  smooth=TRUE)
subsetMainDMRs_ALL = subsetMainDMRs_ALL$table
subsetMainDMRs_ALL$Model = "Subset Main"
subsetMainDMRs_ALL$Region = "ALL"

#####
allRegions = subsetMainDMRs_ALL
save(allRegions, file='rdas/allRegions_DMR_NeuN_Sensitivity.rda')

#####
mergedDMR = rbind(regionalSubset, allRegions )
save(mergedDMR, file='rdas/merged_DMR_NeuN_Sensitivity.rda')
