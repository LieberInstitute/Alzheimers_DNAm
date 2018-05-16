########################
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)

## load data
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')

## get BrNums
BrNums = pd$BrNum
BrNumOriginal = pd$BrNum

pd$BrNum = BrNums

############# 
# fam file ##

### read in fam
# fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10.fam",
    # as.is=TRUE)
fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10.fam",
    as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")
fam$BrNum = fam$FID
fam$BrNum[fam$FID == "Omni2pt5"] = fam$IID[fam$FID == "Omni2pt5"]
# fam[grep("^Br", fam$IID),1:2] = fam[grep("^Br", fam$IID),2:1]
table(pd$BrNum %in% fam$BrNum) # keep samples w/ genotypes
pd$BrNum[!pd$BrNum %in% fam$BrNum]
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M

famOut = fam[which(fam$BrNum %in% pd$BrNum),]

write.table(famOut[,1:2], "genotype_data/samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
setwd('genotype_data')
# bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10"
bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10"
newbfile = "/dcl01/lieber/ajaffe/Steve/Alz/genotype_data/Alz_meth450k_Genotypes_n380_maf005_geno10_hwe1e6"
forceKeep_bfile = "/dcl01/lieber/ajaffe/Steve/Alz/genotype_data/Alz_meth450k_Genotypes_n380_forceKeep"

## extract
system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.1 --maf 0.05 --hwe 0.000001 --make-bed --out", newbfile))

# ## independent and cluster
system(paste("/users/ajaffe/bin/plink --bfile", newbfile, "--indep 100 10 1.25 --out", newbfile))

## MDS components	
system(paste0("/users/ajaffe/bin/plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))	
	
# ## A transpose
system(paste("/users/ajaffe/bin/plink --bfile", newbfile,
	"--recode A-transpose --out", newbfile))

##### Force pull SNPs of interest
system(paste("/users/ajaffe/bin/plink --bfile", bfile, "--keep samples_to_extract.txt --extract range snp_range_to_check.txt --make-bed --out", forceKeep_bfile))
system(paste("/users/ajaffe/bin/plink --bfile", forceKeep_bfile, "--recode A-transpose  --out", forceKeep_bfile))	

################
## read in #####

## read in genotypes
genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")
snp = as.data.frame(genotypes[,-(1:6)])
colnames(snp) = ifelse(grepl("^Br", ss(colnames(snp), "_")),
			ss(colnames(snp), "_"), ss(colnames(snp), "_",2))
						
### update MAP
# load("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10_updatedMap.rda")
load("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10_updatedMap.rda")
snpMap = snpMap[match(genotypes$SNP, snpMap$SNP),]
rownames(snp) = rownames(snpMap) =  snpMap$SNP

## check directionality
swapIndex = which(snpMap$ALT == genotypes$COUNTED)
snp[swapIndex,] = 2-snp[swapIndex,]

### Loading SNPs force kept
forceKeep_geno = read_delim(paste0(forceKeep_bfile, ".traw"), delim="\t")			
forceKeep_snp =	 as.data.frame(forceKeep_geno[,-(1:6)])		
colnames(forceKeep_snp) = ifelse(grepl("^Br", ss(colnames(forceKeep_snp), "_")),
			ss(colnames(forceKeep_snp), "_"), ss(colnames(forceKeep_snp), "_",2))
			
forceKeep_snpMap = read_delim(paste0(forceKeep_bfile, ".bim"), delim="\t", col_names=FALSE)
colnames(forceKeep_snpMap) = c("CHR", "SNP", "CM", "POS", "COUNTED","ALT")
forceKeep_snpMap = as.data.frame(forceKeep_snpMap)
forceKeep_snpMap[,c('Type','newRef','newCount','name','rsNumGuess','chr_hg38','pos_hg38')] <- NA
rownames(forceKeep_snp) = rownames(forceKeep_snpMap) =  forceKeep_snpMap$SNP

#swapIndex = which(forceKeep_snpMap$ALT == forceKeep_geno$COUNTED)
#forceKeep_snp[swapIndex,] = 2-forceKeep_snp[swapIndex,]

pullSNP = forceKeep_snpMap[!rownames(forceKeep_snpMap) %in% rownames(snpMap),,drop=F]

### Merge a single new variant onto the snpMap object
snpMap = rbind(snpMap, forceKeep_snpMap[pullSNP,] )
snp = rbind(snp, forceKeep_snp[pullSNP,] )



#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE)
rownames(mds) = ifelse(grepl("^Br", mds$FID),
			mds$FID, mds$IID)
mds = mds[,-(1:3)]
colnames(mds) = paste0("snpPC",1:ncol(mds))

##########################
## correct BrNumbers #####

## confirm order stayed the same throughout
identical(rownames(mds), colnames(snp))
identical(rownames(mds), pd$BrNum[match(rownames(mds),pd$BrNum)])


#############
## save #####
save(mds, snp, snpMap, compress=TRUE, file = "/dcl01/lieber/ajaffe/Steve/Alz/genotype_data/Alz_meth450k_Genotypes_n380.rda")

##### 
load('/dcl01/lieber/ajaffe/Steve/Alz/genotype_data/Alz_meth450k_Genotypes_n380.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')

pd = cbind(pd, mds[pd$BrNum,])
save(pd,bVals,mVals, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')

### Inspect mds components
pd$genotype_array = famOut[match(pd$BrNum,famOut$FID),'IID']

library(GGally)
library(ggplot2)
pc5_plot = ggpairs(data=as.data.frame(pd[,c('snpPC1','snpPC2','snpPC3','snpPC4','snpPC5','genotype_array') ]) ,
		columns=1:5, # data.frame with variables
        title="First 5 Genotype MDS Components",
		mapping = aes(color = `genotype_array`) )
ggsave(pc5_plot, filename='/dcl01/lieber/ajaffe/Steve/Alz/qc/pairwise_top5_snpPCs_array.pdf',height=40,width=40 )

pc5_plot = ggpairs(data=as.data.frame(pd[,c('snpPC1','snpPC2','snpPC3','snpPC4','snpPC5','Race') ]) ,
		columns=1:5, # data.frame with variables
        title="First 5 Genotype MDS Components",
		mapping = aes(color = `Race`) )
ggsave(pc5_plot, filename='/dcl01/lieber/ajaffe/Steve/Alz/qc/pairwise_top5_snpPCs_race.pdf',height=40,width=40 )

