##
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
library(minfi)
library(limma)
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

dropIndex = which(!pd$keepList)
metadata = pd[-dropIndex,]
metadata$Sample_Plate <- metadata$Sample_Group

metadata[,'Sample name'] = metadata$Sample_Name
metadata[,'title'] = paste0( metadata$Region, "_", metadata$Dx, "_" , metadata$BrNum ) 
metadata[,'source name'] = metadata$Region
metadata[,'organism'] = 'Homo sapiens'
metadata[,'idat file 1'] = paste0(basename(metadata$basePath), "_Grn.idat")
metadata[,'idat file 2'] = paste0(basename(metadata$basePath), "_Red.idat")
#
metadata[,'molecule'] = 'genomic DNA'
metadata[,'label'] = 'Cy3/Cy5'
metadata[,'description'] = "Postmortem human brain tissue samples from four brain regions in an Alzheimers disease case-control series used to generate Human450k DNA methylation data."
metadata[,'platform'] = 'GPL13534'
#
characteristics_vector = c('Sample_Plate','Sentrix_Position','BrNum','Region','Dx','Age','Sex','Race','negControl_PC1','negControl_PC2','NeuN_pos','snpPC1')
colnames(metadata)[colnames(metadata) %in% characteristics_vector] = paste0("characeristics: ", colnames(metadata)[colnames(metadata) %in% characteristics_vector] )

write.csv(metadata[,c('Sample name',
					  'title',
					  'source name',
					  'organism',
					  'idat file 1',
					  'idat file 2',
					  paste0("characeristics: ", characteristics_vector),
					  'molecule',
					  'label',
					  'description',
					  'platform')],'GEO/metaData.csv',row.names=FALSE)

### copy the appropriate idat folders onto a designated directory
cpCall1=paste0(metadata$basePath, "_Red.idat") ### RUN IN BASH
cpCall2=paste0(metadata$basePath, "_Grn.idat") ### RUN IN BASH
cpCall = c(cpCall1, cpCall2)
file.copy(from=cpCall, to=paste0('/dcl01/lieber/ajaffe/Steve/Alz/Paper/GEO/idats'), copy.mode = TRUE, copy.date = FALSE)

### Zip the idat files
zip -r /dcl01/lieber/ajaffe/Steve/Alz/Paper/GEO/idat_files.zip /dcl01/lieber/ajaffe/Steve/Alz/Paper/GEO/idats ### RUN IN BASH

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
bVals <- bVals[,metadata$Sample_Name]

### calculate detection P values
load('rdas/RGset_prefiltered_n394.rda')
RGset = RGset[ ,metadata$Sample_Name]
detection_p.value=detectionP(RGset, type = "m+u")

#load('/home/other/nivanov2/PROJECTS/dura_fibro_DNAm_forPublication/rdas/MSet_normalized_SNPprobes_and_SEXprobes_retained.rda') #MSet_normalized_SNPprobes_and_SEXprobes_retained
#beta=getBeta(MSet_normalized_SNPprobes_and_SEXprobes_retained)

mm=match(rownames(bVals),rownames(detection_p.value))
detection_p.value=detection_p.value[mm,]


all(colnames(detection_p.value)==colnames(bVals)) #TRUE

i=1
p=as.vector(detection_p.value[,i])
mtch=match(colnames(detection_p.value)[i],colnames(bVals))
b=bVals[,mtch]
matrix_processed=data.frame(b,p)
colnames(matrix_processed)=c(paste0(colnames(detection_p.value)[i],'_Beta'),paste0(colnames(detection_p.value)[i],'_Detection_Pval'))

for (i in 2:ncol(detection_p.value)){

	p=detection_p.value[i]
	mtch=match(colnames(detection_p.value)[i],colnames(bVals))

	b=bVals[,mtch]
	tmp=data.frame(b,p)
	colnames(tmp)=c(paste0(colnames(detection_p.value)[i],'_Beta'),paste0(colnames(detection_p.value)[i],'_Detection_Pval'))
	matrix_processed=cbind(matrix_processed,tmp)

}

write.csv(matrix_processed,'GEO/matrix_processed.csv',row.names=TRUE)

###
#detection_p.value=detectionP(RGSet, type = "m+u")

raw=preprocessRaw(RGset) 
unmethRaw=getUnmeth(raw)
methRaw=getMeth(raw)

unmethRaw=unmethRaw[mm,]
methRaw=methRaw[mm,]

all(colnames(unmethRaw)==colnames(methRaw)) #TRUE
all(colnames(detection_p.value)==colnames(methRaw)) #TRUE

all(rownames(detection_p.value)==rownames(unmethRaw)) #TRUE
all(rownames(detection_p.value)==rownames(methRaw)) #TRUE

i=1
p=detection_p.value[,i]
mtch=match(colnames(detection_p.value)[i],colnames(methRaw))
m=methRaw[,mtch]
u=unmethRaw[,mtch]
matrix_signalIntensities=data.frame(u,m,p)
colnames(matrix_signalIntensities)=c(
	paste0(colnames(detection_p.value)[i],'_unmethylated_signal'),
	paste0(colnames(detection_p.value)[i],'_methylated_signal'),
	paste0(colnames(detection_p.value)[i],'_Detection_Pval') )

for (i in 2:ncol(detection_p.value)){

p=detection_p.value[,i]
mtch=match(colnames(detection_p.value)[i],colnames(methRaw))
m=methRaw[,mtch]
u=unmethRaw[,mtch]
temp=data.frame(u,m,p)
colnames(temp)=c(
	paste0(colnames(detection_p.value)[i],'_unmethylated_signal'),
	paste0(colnames(detection_p.value)[i],'_methylated_signal'),
	paste0(colnames(detection_p.value)[i],'_Detection_Pval') )


matrix_signalIntensities=cbind(matrix_signalIntensities,temp)

}

write.csv(matrix_signalIntensities,'GEO/matrix_signal_intensities.csv',row.names=TRUE)

## Subset to samples also present in RNA-seq data
pd_meth = metadata
colnames(pd_meth) <- gsub("characeristics: ", "", colnames(pd_meth))
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_pd_hg38.Rdata')
pd_rna = pd 
pd_rna$Region = plyr::revalue(pd_rna$Region, c("CB"="CRB") )
pd_rna$Dx = plyr::revalue(pd_rna$Dx, c("NC"="Control","AD" = "Alzheimer") )
pd_rna$BrNum = paste0("Br", pd_rna$BRNum)

#### Check for shared brains
pd_meth$matchBrain = paste0(pd_meth$BrNum, "_", pd_meth$Region)
pd_rna$matchBrain = paste0(pd_rna$BrNum, "_", pd_rna$Region)

shared_samples = intersect(pd_meth$matchBrain,pd_rna$matchBrain )
pd_meth = pd_meth[pd_meth$matchBrain%in%shared_samples,]
pd_rna = pd_rna[pd_rna$matchBrain%in%shared_samples,]

### Dropping duplicated brains in RNAseq data
library(dplyr)
dropBrains = pd_rna[duplicated(pd_rna$matchBrain) | duplicated(pd_rna$matchBrain, fromLast=TRUE), ] 
dropBrains = dropBrains[order(dropBrains$BrNum, dropBrains$mitoRate),]
dropBrains = dropBrains[duplicated(dropBrains$matchBrain),'SAMPLE_ID']
pd_rna = pd_rna[!pd_rna$SAMPLE_ID %in% dropBrains, ]

### reorder pd_rna
pd_rna= pd_rna[match(pd_meth$matchBrain,pd_rna$matchBrain), ] 

### Extract count matrix from RNA-sequencing data
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/gene_info.Rdata',verbose=T)

## Filter genes
geneCounts <- geneCounts[geneMap$meanExprs > 0.1, ] # filter with meanExprs >0.1
geneMap <- geneMap[geneMap$meanExprs > 0.1, ]

## Filter people
geneCounts <- geneCounts[,pd_rna$SAMPLE_ID]

#Check then re-name columns of count matrix so they match up w/ RNA-seq data
pd_meth$matchBrain==pd_rna$matchBrain
colnames(geneCounts)<- pd_meth$Sample_Name
## Add a couple columns that are needed to replicate, Gene Symbol, ESMBL gene ID, and gene length
geneCounts = cbind(geneMap[,c('Symbol','Length')],geneCounts)
##
write.csv(geneCounts,file='GEO/geneCounts.csv',row.names=TRUE)

## zip up the text-based files that we'll submit to GEO
cd /dcl01/lieber/ajaffe/Steve/Alz/Paper/GEO ### RUN IN BASH
tar -cvzf GEO_Alzheimers_DNAmData_and_geneCounts_Data.tar.gz idat_files.zip matrix_processed.csv matrix_signal_intensities.csv geneCounts.csv readme.txt GA_illumina_methylation_AlzheimersDisease_FourRegions_DNAm.xls ### RUN IN BASH
