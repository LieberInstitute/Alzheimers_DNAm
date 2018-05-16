#
library(minfi)
setwd('/dcl01/lieber/ajaffe/Steve/Alz/')
##
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/MSet_SQN_postfiltered.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/processed_data_postfiltered.rda')

##
dropList1=openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/qc/n398_meth450k_genotype_problematic_samples.xlsx',sheet=1)
dropList1$Decision='Drop'
dropList2=openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/qc/n398_meth450k_genotype_problematic_samples.xlsx',sheet=2)

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

#length(keep_genotyped)
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

### Save out final sample set used for meQTL analysis
save(pd,bVals,mVals,file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
save(pd,Mset_SQN,file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_Mset_SQN_postfiltered.rda')