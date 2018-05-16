#qsub -l mf=40G,h_vmem=80G,h_stack=256M,h_fsize=50G -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 00_load_idats.R

library(minfi)
library(jaffelab)

#00 loading in methylation data 
###############################
sampleSheet = read.csv("/dcl01/lieber/ajaffe/NEXSAN2/ajaffe/450k/BrainDataGSK/samplePlateInfo_GSK.csv",
	as.is=TRUE)
pheno = read.csv("/users/ajaffe/Lieber/Projects/450k/alz/samples_for_gsk_450k.csv", as.is=TRUE)
colnames(pheno)[c(1,4)] = c("BrNum","Sample_Name")

pd = merge(sampleSheet, pheno, by="Sample_Name",all.x=TRUE)

### only brain
pd = pd[pd$Region %in% c("cerebellum", "DLPFC", "ERC", "Hippo"),]
pd$BrNum = paste0("Br", pd$BrNum)

## case control status
pheno2 = read.csv("/users/ajaffe/Lieber/Projects/450k/alz/gsk_brain_phenos.csv", as.is=TRUE)

m1 = match(pd$BrNum, paste0("Br",pheno2$BrNumOld))
m2 = match(pd$BrNum, pheno2$BrNum1)
m1[is.na(m1)] = m2[is.na(m1)]

pd = cbind(pd, pheno2[m1,6:9])
pd = pd[!is.na(pd$Dx),]

## format
pd$Region[grep("cerebellum",pd$Region)] = "CRB"
pd$Region = toupper(pd$Region)

pd$basePath = paste0("/dcl01/lieber/ajaffe/NEXSAN2/ajaffe/450k/BrainDataGSK/",
	pd$Sentrix_ID, "/", pd$Sentrix_ID,"_", pd$Sentrix_Position)
	
#### update brain numbers
pheno3 = read.csv("/users/ajaffe/Lieber/Projects/450k/alz/pheno_info_from_straub.csv", as.is=TRUE,row.names=1)
mm = match(pd$BrNum, paste0("Br", pheno3$OLD_BrNUM))
pd$BrNum[!is.na(mm)] = pheno3$BrNum[mm[!is.na(mm)]]
pd$BrNum = gsub("Br","", pd$BrNum)

## drop those not on list
pd = pd[pd$BrNum %in% pheno3$BrNum,]

#### read in the data
RGset = read.metharray(pd$basePath)
pData(RGset) <- DataFrame(pd)
save(RGset, pd, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/RGset_n398.rda')