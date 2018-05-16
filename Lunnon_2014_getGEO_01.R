library(minfi)
library(GEOquery)
#GSE43414

## Download raw data
#download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59685/soft/',destfile='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_SOFT')

##
lunnon = getGEO('GSE59685')
lunnon = lunnon[[1]]

processed_bVals <- as.data.frame(exprs(lunnon))
pd = pData(phenoData(lunnon))
pd= pd[,c('title',"geo_accession", grep("characteristics",colnames(pd),value=TRUE), grep("data_processing",colnames(pd),value=TRUE),'relation') ]

pd = plyr::rename(pd, c("characteristics_ch1" = "SampleID", "characteristics_ch1.1"="Chip", "characteristics_ch1.2" ="Dx", "characteristics_ch1.3" = "braak_stage", "characteristics_ch1.4" = "Sex", "characteristics_ch1.5" = "Blood Age", "characteristics_ch1.6"="Age","characteristics_ch1.7" ="Region") )
rownames(pd) = NULL
pd=pd[,-grep("data_processing",colnames(pd))]

pd$SampleID = gsub("subjectid: ", "", pd$SampleID)

pd$Chip = gsub("barcode: ", "", pd$Chip)

pd$Dx = gsub("ad.disease.status: ", "", pd$Dx)
pd$Dx = plyr::revalue(pd$Dx, c("C"="Control","Exclude"=NA,"AD"="Alzheimers") )
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimers") )

pd$braak_stage = gsub("braak.stage: ", "", pd$braak_stage)
pd$braak_stage = plyr::revalue(pd$braak_stage, c("Exclude"=NA) )
pd$braak_stage = as.numeric(pd$braak_stage)

pd$Sex = gsub("Sex: ", "", pd$Sex)
pd$Sex = plyr::revalue(pd$Sex, c("FEMALE"="F", "MALE"="M") )
pd$Sex = factor(pd$Sex, levels= c("M","F") )

pd$Blood_Age = gsub("age.blood: ", "", pd$Blood_Age)
pd$Blood_Age[pd$Blood_Age=="NA"] <- NA
pd$Blood_Age = as.numeric(pd$Blood_Age)

pd$Age = gsub("age.brain: ", "", pd$Age)
pd$Age[pd$Age=="NA"] <- NA
pd$Age = as.numeric(pd$Age)

pd$Region = gsub("source tissue: ", "", pd$Region)
pd$Region = plyr::revalue(pd$Region, c("cerebellum"="CRB", "superior temporal gyrus" = "STG", "whole blood"= "blood","entorhinal cortex" = 'ERC',  "frontal cortex" ="PFC" ) )

save(pd,file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/raw_signalIntensity.rda')

##
getGEOSuppFiles('GSE59685', makeDirectory = TRUE, baseDir = '/dcl01/lieber/ajaffe/Steve/Alz/')
sapply(list.files("/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/", pattern= "gz", full.names=T), gunzip)
untar("/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_RAW.tar", exdir = "/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/raw")

##
bVals = data.table::fread('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_betas.csv') 
processed_bVals=bVals[-(1:2),]
processed_bVals = as.data.frame(processed_bVals)
rownames(processed_bVals) <- processed_bVals[,1]
processed_bVals = processed_bVals[,-1]
processed_bVals= processed_bVals[,pd$Chip]

rn=rownames(processed_bVals)
processed_bVals = sapply(processed_bVals, as.numeric)
rownames(processed_bVals)=rn

###
signalIntensities = data.table::fread('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_signal_intensities.csv') 
write.csv(signalIntensities, file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_signal_intensities_preprocessed.csv',row.names=F)

lumiMethyl <- methylumi::methylumiR('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_signal_intensities_preprocessed.csv',sep=',')

MSet = minfi::readGEORawFile(file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_signal_intensities_preprocessed.csv',Uname='Unmethylated Signal', Mname='Methylated Signal')

#just annotation
#GPL13534 = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GPL13534_HumanMethylation450_15017482_v.1.1.csv')


####
save(pd,signalIntensities, file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/raw_signalIntensity.rda')
save(pd,processed_bVals, file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/processed_bVals.rda')
save(lunnon, file='/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/gse_object.rda')