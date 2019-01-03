setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
## Check for probe-QC vs. AD GWAS association
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')
###

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

###
bVals_Without_Remapping_Drop <- bVals
bVals <- bVals[-drop_hg38_unmappable, ] 



################## Checking loss of probes in GWAS loci ####################
### annotation
#load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation_subset.rda',verbose=T) #load hg38 position annotation
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

######### USING LD BLOCKS INSTEAD OF 20KB FLANKER
ad_ld_friends <- read.csv('/dcl01/lieber/ajaffe/Steve/Alz/raggR_21SNP_AD_Meta_LD_friends_r2_60.csv')
ad_ld_friends = ad_ld_friends[!duplicated(ad_ld_friends$SNP2.Name),]
ad_ld_friends$SNP1.Chr = paste0('chr', ad_ld_friends$SNP1.Chr )
ad_ld_friends$SNP2.Chr = paste0('chr', ad_ld_friends$SNP2.Chr ) 

############ build the risk loci #######
library(GenomicRanges)
risk_loci_grl = makeGRangesFromDataFrame(ad_ld_friends, start.field='SNP2.Pos',end.field='SNP2.Pos',seqnames.field='SNP2.Chr', keep.extra.columns=T )
risk_loci_grl = split(risk_loci_grl, as.factor(risk_loci_grl$SNP1.Name))
risk_loci_grl = lapply(risk_loci_grl, range )
risk_loci_grl = GRangesList(risk_loci_grl)
#risk_loci_gr =  do.call(getMethod(c, "GenomicRanges"),risk_loci_grl)
#do.call("c",risk_loci_grl)
risk_loci_gr = unlist(risk_loci_grl, recursive = TRUE, use.names = TRUE)
#risk_loci_gr=unlist(risk_loci_grl)
names(risk_loci_gr) <- names(risk_loci_grl)
risk_loci_gr$coordinates=paste0(as.character(seqnames(risk_loci_gr)),":",start(risk_loci_gr),"-",end(risk_loci_gr) )
ld_loci_coverage = sum(width(risk_loci_gr ) )


########## overlay risk regions with cpgs present #########
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
risk_loci_list = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_loci.csv')
goldsetGR=GRanges(ranges=IRanges(start=ann450k$pos, end=ann450k$pos), seqnames=ann450k$chr)
names(goldsetGR) <- ann450k$Name

overlaps = findOverlaps(risk_loci_gr, goldsetGR)

ann450k$Overlaps_AD_GeneticRiskLoci <- FALSE
ann450k$Overlaps_AD_GeneticRiskLoci[subjectHits(overlaps)] <- TRUE

ann450k[,c('SNP','loci_coordinates') ] <- NA
ann450k[subjectHits(overlaps),'SNP'] = names(risk_loci_gr)[queryHits(overlaps)]
ann450k[subjectHits(overlaps),'loci_coordinates'] =  risk_loci_gr[queryHits(overlaps)]$coordinates

ad_genetic_risk_loci = ann450k[ann450k$Overlaps_AD_GeneticRiskLoci, c('Name','SNP','loci_coordinates')]
ad_genetic_risk_loci$indexSnp_rsNum <- jaffelab::ss(ad_genetic_risk_loci$SNP,":",1)
##
kk = match(ad_genetic_risk_loci$indexSnp_rsNum, risk_loci_list$SNP)
ad_genetic_risk_loci$indexSnp_Chr_hg19 = risk_loci_list$Chr[kk]
ad_genetic_risk_loci$indexSnp_Pos_hg19 = risk_loci_list$Position[kk]
ad_genetic_risk_loci$indexSnp_nearestGene = risk_loci_list$Closest.gene[kk]
dim(ad_genetic_risk_loci)

dropped_gwas_probes = sum(!rownames(ad_genetic_risk_loci)%in%rownames(bVals) )
tested_gwas_probes = sum(rownames(ad_genetic_risk_loci)%in%rownames(bVals) )
tested_notGwas_probes = nrow(bVals)-tested_gwas_probes
dropped_notGwas_probes = nrow(ann450k) - dropped_gwas_probes - tested_gwas_probes - tested_notGwas_probes

tab = matrix(c(dropped_gwas_probes,tested_gwas_probes,dropped_notGwas_probes,tested_notGwas_probes ),nrow=2 )
rownames(tab) = c('Not tested', 'Tested')
colnames(tab) = c('Within GWAS loci', 'Outside GWAS loci')
tab
fisher.test(tab)
table(as.data.frame(ad_genetic_risk_loci[!rownames(ad_genetic_risk_loci)%in%rownames(bVals), ] )$indexSnp_nearestGene)
as.data.frame(ad_genetic_risk_loci[!rownames(ad_genetic_risk_loci)%in%rownames(bVals), ] )


### hg38 remapping probes
hg38_remap_drop = ann450k[!rownames(ann450k)%in%rownames(goldset),'Name']

hg38_NotDropped_gwas_probes = sum(!rownames(ad_genetic_risk_loci)%in%hg38_remap_drop )
hg38_Dropped_gwas_probes = sum(rownames(ad_genetic_risk_loci)%in%hg38_remap_drop )
hg38_Dropped_NotGwas_probes = sum(!rownames(goldset)%in%rownames(ad_genetic_risk_loci) )
hg38_NotDropped_NotGwas_probes =  nrow(ann450k) - hg38_NotDropped_gwas_probes - hg38_Dropped_gwas_probes -  hg38_Dropped_NotGwas_probes

tab = matrix(c(hg38_Dropped_gwas_probes,hg38_NotDropped_gwas_probes,hg38_NotDropped_NotGwas_probes,hg38_Dropped_NotGwas_probes ),nrow=2 )
rownames(tab) = c('Hg38 Dropped', 'Not hg38 Dropped')
colnames(tab) = c('Within GWAS loci', 'Outside GWAS loci')
tab
fisher.test(tab)
as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),hg38_remap_drop), ] )
table(as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),hg38_remap_drop), ] )$indexSnp_nearestGene)



### xReactiveProbes
xReactiveProbes <- read.csv(file="/users/ssemick/Tox_Expr/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
drop_xreactive <- xReactiveProbes$TargetID

xReactive_NotDropped_gwas_probes = sum(!rownames(ad_genetic_risk_loci) %in% drop_xreactive )
xReactive_Dropped_gwas_probes = sum(rownames(ad_genetic_risk_loci)%in%drop_xreactive )
xReactive_Dropped_NotGwas_probes = length(drop_xreactive) - xReactive_Dropped_gwas_probes
xReactive_NotDropped_NotGwas_probes =  nrow(ann450k) - xReactive_NotDropped_gwas_probes - xReactive_Dropped_gwas_probes -  xReactive_Dropped_NotGwas_probes

tab = matrix(c(xReactive_Dropped_gwas_probes,xReactive_NotDropped_gwas_probes,xReactive_Dropped_NotGwas_probes,xReactive_NotDropped_NotGwas_probes ),nrow=2 )
rownames(tab) = c('xReative Dropped', 'Not xReactive Dropped')
colnames(tab) = c('Within GWAS loci', 'Outside GWAS loci')
tab
fisher.test(tab)
as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),drop_xreactive), ] )
table(as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),drop_xreactive), ] )$indexSnp_nearestGene)


### SBE and CpG probes
load('rdas/RGset_prefiltered_n394.rda')
library(minfi)

SNP_info = getSnpInfo(RGset, snpAnno = NULL)
Drop_SNP = rownames (SNP_info[which( (SNP_info$SBE_maf >= 0.01 | SNP_info$CpG_maf >= 0.01) ),] )

SNP_NotDropped_gwas_probes = sum(!rownames(ad_genetic_risk_loci)%in%Drop_SNP )
SNP_Dropped_gwas_probes = sum(rownames(ad_genetic_risk_loci)%in%Drop_SNP )
SNP_Dropped_NotGwas_probes = length(Drop_SNP) - SNP_Dropped_gwas_probes
SNP_NotDropped_NotGwas_probes =  nrow(ann450k) - SNP_NotDropped_gwas_probes - SNP_Dropped_gwas_probes -  SNP_Dropped_NotGwas_probes

tab = matrix(c(SNP_Dropped_gwas_probes,SNP_NotDropped_gwas_probes,SNP_Dropped_NotGwas_probes,SNP_NotDropped_NotGwas_probes ),nrow=2 )
rownames(tab) = c('SNP Dropped', 'Not SNP Dropped')
colnames(tab) = c('Within GWAS loci', 'Outside GWAS loci')
tab
fisher.test(tab)
as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),Drop_SNP), ] )
table(as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),Drop_SNP), ] )$indexSnp_nearestGene)


### Quality control
#
detP <- detectionP(RGset)
failed <- detP > 0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
keep <- (rowMeans(failed)<0.05)  #removes probes that fail in over 5 percent of samples
table(keep) #5% rule "An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform"
#
QC_Drop = names(keep)[!keep]

QC_NotDropped_gwas_probes = sum(!rownames(ad_genetic_risk_loci)%in%QC_Drop )
QC_Dropped_gwas_probes = sum(rownames(ad_genetic_risk_loci)%in%QC_Drop )
QC_Dropped_NotGwas_probes = length(QC_Drop) - QC_Dropped_gwas_probes
QC_NotDropped_NotGwas_probes =  nrow(ann450k) - QC_NotDropped_gwas_probes - QC_Dropped_gwas_probes -  QC_Dropped_NotGwas_probes

tab = matrix(c(QC_Dropped_gwas_probes,QC_NotDropped_gwas_probes,QC_Dropped_NotGwas_probes,QC_NotDropped_NotGwas_probes ),nrow=2 )
rownames(tab) = c('QC Dropped', 'Not QC Dropped')
colnames(tab) = c('Within GWAS loci', 'Outside GWAS loci')
tab
fisher.test(tab)

as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),QC_Drop), ] )
table(as.data.frame(ad_genetic_risk_loci[intersect(rownames(ad_genetic_risk_loci),QC_Drop), ] )$indexSnp_nearestGene)
