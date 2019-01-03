#
library(minfi)
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
##
load('rdas/MSet_SQN_postfiltered.rda')
load('rdas/processed_data_postfiltered.rda')

##
dropList1=openxlsx::read.xlsx('qc/n398_meth450k_genotype_problematic_samples.xlsx',sheet=1)
dropList1$Decision='Drop'
dropList2=openxlsx::read.xlsx('qc/n398_meth450k_genotype_problematic_samples.xlsx',sheet=2)

## Drop samples which matched the wrong genotype
dropList1 = dropList1[dropList1$Decision=="Drop",'MethSampleName' ]
dropId1=which(pd$Sample_Name %in% dropList1)
length(dropId1) #5
length(unique(pd[dropId1,'BrNum'])) #2
#Br5083 and Br1504
#

pd = pd[-dropId1,]
bVals = bVals[,-dropId1]
mVals = mVals[,-dropId1]
Mset_SQN = Mset_SQN[, -dropId1]

## Drop samples which did not match their own genotype
dropList2 = dropList2$Sample_Name
dropId2=which(pd$Sample_Name %in% dropList2)
length(dropId2) #4
length(unique(pd[dropId2,'BrNum'])) #Br1060 x4
#
pd = pd[-dropId2,]
bVals = bVals[,-dropId2]
mVals = mVals[,-dropId2]
Mset_SQN = Mset_SQN[, -dropId2]

dim(pd) #380
## Check
table(pd$Region,pd$Dx)

## Drop cerebellum outliers
dropId3 = which(pd$Region=="CRB" & pd$BrNum%in%c("Br1909","Br1615","Br2257"))
 
pd = pd[-dropId3,]
bVals = bVals[,-dropId3]
mVals = mVals[,-dropId3]
Mset_SQN = Mset_SQN[, -dropId3]

dim(pd) #377

### Check colnames
all(pd$Sample_Name ==colnames(bVals))
all(pd$Sample_Name==colnames(mVals))
all(pd$Sample_Name==colnames(Mset_SQN))

### Add info about keep alz keep list
keep_list = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/keep_list.csv')
keep_list[!is.na(keep_list)] <- paste0("Br", as.character(as.numeric(keep_list[!is.na(keep_list)] )))

pd$keepList = ifelse(pd$Dx=="Control",NA,FALSE)
pd[pd$Region=="CRB" & pd$BrNum%in%keep_list$CB,'keepList' ] <- TRUE
pd[pd$Region=="DLPFC" & pd$BrNum%in%keep_list$DLPFC,'keepList' ] <- TRUE
pd[pd$Region=="ERC" & pd$BrNum%in%keep_list$ERC,'keepList' ] <- TRUE
pd[pd$Region=="HIPPO" & pd$BrNum%in%keep_list$Hippo,'keepList' ] <- TRUE

all_regions = table(pd$BrNum)==4
all_regions=names(all_regions)[all_regions]
pd$allRegions = pd$BrNum %in% all_regions

pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd$keepList[is.na(pd$keepList)] = TRUE
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'

### Save out final sample set used for meQTL analysis
#Load in pheno data from Rahul 9/12/2018
updated_keep_list = openxlsx::read.xlsx('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis/NIMH and LIBD IDs for AD CASES.xlsx',start=2)
updated_keep_list$BrNumerical = as.numeric( gsub("BR","",updated_keep_list$LBR.Brain.ID) )
###
dat = cbind(pd, updated_keep_list[match( as.numeric (gsub("Br","",pd$BrNum )), updated_keep_list$BrNumerical  ),] ) 
oldbrnums = read.csv('/users/ajaffe/Lieber/Projects/450k/alz/pheno_info_from_straub.csv')
dat = cbind(dat, Old_BrNum = oldbrnums[match( as.numeric (gsub("Br","",pd$BrNum )), oldbrnums$BrNum),'OLD_BrNUM'] )
table(dat$DxOrdinal, dat$AD.RATING)

##
info_out = dat[!duplicated(dat[,c('BrNum','DxOrdinal')]),c('BrNum','Old_BrNum','Brain.ID','AD.RATING','Dx','DxOrdinal')]
info_out$INCORRECTLY_USED_IN_ORDINAL_MODEL = ifelse(!is.na(info_out$`AD.RATING`) & info_out$DxOrdinal=='Alz Drop', 'Yes','No')
write.csv(info_out, file='/dcl01/lieber/ajaffe/Steve/Alz/checking_ordinal_model.csv', row.names=F)

BrNums= dat[!(as.numeric (gsub("Br","",dat$BrNum ) )%in%updated_keep_list$BrNumerical | dat$DxOrdinal=="Control"),c('BrNum','Brain.ID')]
BrNums=data.frame(BrNum=unique(BrNums))
write.csv(BrNums,file='/dcl01/lieber/ajaffe/Steve/Alz/BrNums_for_Intermediate_Group.csv')
###
save(pd,bVals,mVals,file='rdas/cleanSamples_n377_processed_data_postfiltered.rda')
save(pd,Mset_SQN,file='rdas/cleanSamples_n377_Mset_SQN_postfiltered.rda')