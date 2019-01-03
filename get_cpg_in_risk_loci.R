#Get CpGs in genetic risk loci
risk_loci_list = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_loci.csv')

############ build the risk loci #######
library(GenomicRanges)
risk_loci_list$start = risk_loci_list$Position
risk_loci_list$end = risk_loci_list$Position
risk_loci_gr = GRanges(risk_loci_list)
risk_loci_gr = reduce(flank(risk_loci_gr,both=TRUE,width=5e4) ) #50kb flank on each side of genetic marker

########## overlay risk regions #########
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
goldsetGR=GRanges(ranges=IRanges(start=goldset$pos, end=goldset$pos), seqnames=goldset$chr)
names(goldsetGR) <- goldset$Name

overlaps = findOverlaps(risk_loci_gr, goldsetGR)

goldset$Overlaps_AD_GeneticRiskLoci <- FALSE
goldset$Overlaps_AD_GeneticRiskLoci[subjectHits(overlaps)] <- TRUE

goldset[,c('snpRsNum','snpChr_hg19','snpPos_hg19','Stage1_MetaP','Stage2_MetaP','Overall_MetaP','nearestGeneToSNP') ] <- NA
goldset[subjectHits(overlaps), c('snpRsNum','snpChr_hg19','snpPos_hg19','Stage1_MetaP','Stage2_MetaP','Overall_MetaP','nearestGeneToSNP')] = risk_loci_list[queryHits(overlaps), c('SNP','Chr','Position', 'Stage1_MetaP','Stage2_MetaP','Overall_MetaP','Closest gene') ]

risk_loci_list[queryHits(overlaps), c('Position') ] - goldset[subjectHits(overlaps), c('pos')]

goldset$SNP_Distance_from_CpG_hg19 =  goldset$pos - goldset$snpPos_hg19

ad_genetic_risk_loci = goldset[goldset$Overlaps_AD_GeneticRiskLoci, c('Name','snpRsNum','snpChr_hg19','snpPos_hg19','Stage1_MetaP','Stage2_MetaP','Overall_MetaP','nearestGeneToSNP')]

write.csv(ad_genetic_risk_loci,file='/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_loci_cpgs.csv',row.names=F)
flank_loci_coverage = sum(width(risk_loci_gr))
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
risk_loci_gr =  unlist(as(risk_loci_grl, "GRangesList") ) #do.call(getMethod(c, "GenomicRanges"),risk_loci_grl)
names(risk_loci_gr) <- names(risk_loci_grl)
risk_loci_gr$coordinates=paste0(as.character(seqnames(risk_loci_gr)),":",start(risk_loci_gr),"-",end(risk_loci_gr) )

ld_loci_coverage = sum(width(risk_loci_gr ) )

########## overlay risk regions with cpgs present #########
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
goldsetGR=GRanges(ranges=IRanges(start=goldset$pos, end=goldset$pos), seqnames=goldset$chr)
names(goldsetGR) <- goldset$Name

overlaps = findOverlaps(risk_loci_gr, goldsetGR)

goldset$Overlaps_AD_GeneticRiskLoci <- FALSE
goldset$Overlaps_AD_GeneticRiskLoci[subjectHits(overlaps)] <- TRUE

goldset[,c('SNP','loci_coordinates') ] <- NA
goldset[subjectHits(overlaps),'SNP'] = names(risk_loci_gr)[queryHits(overlaps)]
goldset[subjectHits(overlaps),'loci_coordinates'] =  risk_loci_gr[queryHits(overlaps)]$coordinates

ad_genetic_risk_loci = goldset[goldset$Overlaps_AD_GeneticRiskLoci, c('Name','SNP','loci_coordinates')]
ad_genetic_risk_loci$indexSnp_rsNum <- jaffelab::ss(ad_genetic_risk_loci$SNP,":",1)
##
kk = match(ad_genetic_risk_loci$indexSnp_rsNum, risk_loci_list$SNP)
ad_genetic_risk_loci$indexSnp_Chr_hg19 = risk_loci_list$Chr[kk]
ad_genetic_risk_loci$indexSnp_Pos_hg19 = risk_loci_list$Position[kk]
ad_genetic_risk_loci$indexSnp_nearestGene = risk_loci_list$Closest.gene[kk]
dim(ad_genetic_risk_loci)

write.csv(ad_genetic_risk_loci,file='/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_LD_block_loci_cpgs.csv',row.names=F)
