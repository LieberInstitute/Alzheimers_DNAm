# Script to download public DNAm data from: DNA methylation analysis on purified neurons and glia dissects age and Alzheimer’s disease-specific changes in the human cortex

#load libraries
library(minfi)
library(GEOquery)
 
#pull processe GEO files + phenotype info
gasparoni = getGEO('GSE66351')
gasparoni = gasparoni[[1]]

#pull raw idats
#getGEOSuppFiles('GSE66351', makeDirectory = TRUE, baseDir = '/dcl01/lieber/ajaffe/Steve/Alz/')
setwd('/dcl01/lieber/ajaffe/Steve/Alz/GSE66351')

## get phenotype info
pd = pData(phenoData(gasparoni))
pd= pd[,c('title',"geo_accession", grep("characteristics",colnames(pd),value=TRUE)) ]

## clean pd info
pd = plyr::rename(pd, c("characteristics_ch1" = "cell_type", "characteristics_ch1.1"="Dx", "characteristics_ch1.2" = "braak_stage", "characteristics_ch1.3" = "Region", "characteristics_ch1.4" = "Age", "characteristics_ch1.5"="Sex","characteristics_ch1.6" ="donor_id","characteristics_ch1.7"="Sentrix_ID","characteristics_ch1.8"="Sentrix_Position") )

rownames(pd) = NULL

pd$cell_type = gsub("cell type: ", "", pd$cell_type)
pd$cell_type = factor(pd$cell_type, levels=c("bulk","Glia","Neuron") )

pd$Dx = gsub("diagnosis: ", "", pd$Dx)
pd$Dx = plyr::revalue(pd$Dx, c("CTRL"="Control","AD"="Alzheimers") )
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimers") )

pd$braak_stage = gsub("braak_stage: ", "", pd$braak_stage)
pd$braak_stage = as.numeric(pd$braak_stage)

pd$Region = gsub("brain_region: ", "", pd$Region)
pd$Region = gsub("brain_region: ", "", pd$Region)

pd$Age = gsub("age: ", "", pd$Age)
pd$Age = as.numeric(pd$Age)

pd$Sex = gsub("Sex: ", "", pd$Sex)
pd$Sex = factor(pd$Sex, levels= c("M","F") )

pd$donor_id = gsub("donor_id: ", "", pd$donor_id)

pd$Sentrix_ID = gsub("sentrix_id: ", "", pd$Sentrix_ID)
pd$Sentrix_Position = gsub("sentrix_position: ", "", pd$Sentrix_Position)
###
pd$basePath = paste0("/dcl01/lieber/ajaffe/Steve/Alz/GSE66351/idats/",
	pd$geo_accession, "_", pd$Sentrix_ID,"_", pd$Sentrix_Position)
	
RGset = read.metharray(pd$basePath)
pData(RGset) <- DataFrame(pd)
save(RGset, pd, file='rdas/RGset_n190.rda')