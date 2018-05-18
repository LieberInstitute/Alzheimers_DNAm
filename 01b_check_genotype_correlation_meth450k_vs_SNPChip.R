library(jaffelab)
library(VariantAnnotation)
library(LIBDpheno)
library(pheatmap)
library(RColorBrewer)
library(minfi)
col.pal = brewer.pal(9,"Blues")
###########################
### oberved genotypes #####

## those genotyped SNPs
genotyped = readVcf("/dcl01/lieber/ajaffe/Steve/BrainSwap_Methylation/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_GenotypingBarcode_Methylation450k.vcf.gz", genome="hg19")

## get sex info
fam = read.table("/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10.fam",
	as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")

## fix 5 digit Brain numbers
fam$newBrNum = fam$FID
fam$newBrNum[nchar(fam$newBrNum)==7] = paste0("Br", 
	substr(fam$newBrNum[nchar(fam$newBrNum)==7], 4, 7))
fam$SampleID = paste0(fam$FID, "_", fam$IID)
colnames(genotyped)  = fam$SampleID
colData(genotyped)$BrNum = ss(colnames(genotyped), "_")

######
snpsGeno = geno(genotyped)$GT
snpsGeno[snpsGeno == "."] = NA
snpsGeno[snpsGeno == "0/0"] = 0
snpsGeno[snpsGeno == "0/1"] = 1
snpsGeno[snpsGeno == "1/1"] = 2
class(snpsGeno) = "numeric"

############ methylation called genotypes ##########
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/RGset_n398.rda')
pd$BrNum[!is.na(pd$BrNum)]= paste0("Br", as.character(as.numeric(gsub("Br", "", pd$BrNum[!is.na(pd$BrNum)]) )) )
pd$Chip = paste0(pd$Sentrix_ID,"_",pd$Sentrix_Position)
#Change Hippocampus to Hippo
sampleNames(RGset) <- pd$Chip
snps <- getSnpBeta(RGset)
colnames(snps) = paste0(pd[match(colnames(snps),pd$Chip),'BrNum'],"_",pd[match(colnames(snps),pd$Chip),'Sample_Name'])
snps = snps[ ,!is.na(colnames(snps))]

#Subset to shared SNPs
shared_snps = intersect( rownames(snps), rownames(snpsGeno ) )
snps=snps[shared_snps,]
snpsGeno=snpsGeno[shared_snps,]

## Call genotypes from beta clouds (k-means algorithm )
calculate.beta.genotypes <- function(snp.betas, centers=c(0.15,0.5,0.85)) {
    x <- t(apply(snp.betas,1,function(x) {
        tryCatch(kmeans(x, centers=centers)$cluster - 1,
                 error=function(e) {
                     cluster <- rep(1,ncol(snp.betas))
                     cluster[which(x < min(centers))] <- 0
                     cluster[which(x > max(centers))] <- 2
                     cluster
                 })
    }))
    dimnames(x) <- dimnames(snp.betas)
    x
}
beta.genotypes <- calculate.beta.genotypes(snps)

### Subset to group sharing genotypes and methylation to find which SNPs to flip
shared_people = intersect(jaffelab::ss(colnames(beta.genotypes),"_"), jaffelab::ss(colnames(snpsGeno),"_"))

beta.genotypes.sharedSubset <- beta.genotypes[,jaffelab::ss(colnames(beta.genotypes),"_")%in%shared_people]
colnames(beta.genotypes.sharedSubset) <- jaffelab::ss(colnames(beta.genotypes.sharedSubset),"_")
beta.genotypes.sharedSubset = beta.genotypes.sharedSubset[,!duplicated(colnames(beta.genotypes.sharedSubset))]

snpsGeno.sharedSubset <- snpsGeno[,jaffelab::ss(colnames(snpsGeno),"_")%in%shared_people]
snpsGeno.sharedSubset = snpsGeno.sharedSubset[,match(colnames(beta.genotypes.sharedSubset),jaffelab::ss(colnames(snpsGeno.sharedSubset),"_") )] # fix col order

colnames(snpsGeno.sharedSubset) <- jaffelab::ss(colnames(snpsGeno.sharedSubset),"_")
snpsGeno.sharedSubset <- snpsGeno.sharedSubset[, !duplicated(jaffelab::ss(colnames(snpsGeno.sharedSubset),"_"))] #drop unnecessary duplicates, need same matrix dimensions
colnames(beta.genotypes.sharedSubset) == jaffelab::ss(colnames(snpsGeno.sharedSubset),"_")#checking

### Find SNPs to flip (0=2, 2=0)
    snp.dosageFlip <- sapply(shared_snps, function(snp) {
        beta.factor <- factor(beta.genotypes.sharedSubset[snp,],levels=0:2)
        geno.factor <- factor(snpsGeno.sharedSubset[snp,],levels=0:2)
        counts <- table(beta.factor, geno.factor)
        diag1 <- sum(counts[1,1] + counts[2,2] + counts[3,3])
        diag2 <- sum(counts[3,1] + counts[2,2] + counts[1,3])
        if (diag1 < diag2) {
			TRUE
        }
        else
            FALSE
    })
snp.dosageFlip=names(snp.dosageFlip)[snp.dosageFlip]
	
## Flip the SNPs for methylation identifed above
beta.genotypes.flipped = beta.genotypes
beta.genotypes.flipped[snp.dosageFlip, ] = 2 - beta.genotypes.flipped[snp.dosageFlip, ]
snps.flipped = snps
snps.flipped[snp.dosageFlip, ] = -1*snps.flipped[snp.dosageFlip, ]

### Checking that flipped genos make more sense
beta.genotypes.flipped_subset = beta.genotypes.flipped[,jaffelab::ss(colnames(beta.genotypes.flipped),"_")%in%shared_people]
colnames(beta.genotypes.flipped_subset) <- jaffelab::ss(colnames(beta.genotypes.flipped_subset),"_")
beta.genotypes.flipped_subset = beta.genotypes.flipped_subset[,!duplicated(colnames(beta.genotypes.flipped_subset))]

colnames(beta.genotypes.flipped_subset) == jaffelab::ss(colnames(snpsGeno.sharedSubset),"_")#checking
rowSums(abs(beta.genotypes.flipped_subset-snpsGeno.sharedSubset),na.rm=TRUE)
table(rowSums(abs(beta.genotypes.flipped_subset-snpsGeno.sharedSubset),na.rm=TRUE)>100)

## Calculate correlation now that alleles have been correctly flipped
snpCor = cor(snpsGeno,beta.genotypes.flipped, use="pairwise.comp")
snpCor2 = cor(snpsGeno,snps.flipped, use="pairwise.comp")
## Correlation across shared samples
corShared = cor(beta.genotypes.flipped_subset,snpsGeno.sharedSubset, use="pairwise.comp")


summary(diag(abs(corShared[jaffelab::ss(rownames(corShared),"_") %in% shared_people, ])))

### Samples with high correlation, unmatching BrNums
missingGenotypes = pd[!pd$BrNum%in%jaffelab::ss(rownames(snpCor),"_",1),c('BrNum')] #these samples do not have genotypes at all
 #write.table(missingGenotypes,file='BrNum_without_genotype_data.txt', sep="\n",row.names=FALSE,col.names=FALSE)
snpCorMat=snpCor
matchInd = as.data.frame(which(snpCorMat > 0.7, arr.ind=TRUE, useNames=TRUE))
matchInd$GenoInfo = rownames(snpCorMat)[matchInd$row]
matchInd$GenoBrNum = ss(matchInd$GenoInfo, "_")
matchInd$MethInfo = colnames(snpCorMat)[matchInd$col]
matchInd$MethBrNum = ss(matchInd$MethInfo, "_")
matchInd$MethSampleName =  paste0(ss(matchInd$MethInfo, "_",2) )
matchInd$MethRegion = pd[match( matchInd$MethSampleName, pd$Sample_Name),'Region']

## add correlations
for(i in 1:nrow(matchInd)) matchInd$genoCor[i] = snpCorMat[matchInd$row[i], matchInd$col[i]]
matchInd = matchInd[,c(3, 4,5,6,7,8,9)]
rownames(matchInd) <- NULL

# check duplicates
mismatched_geno = matchInd[matchInd$MethInfo %in% matchInd$MethInfo[
	matchInd$MethBrNum != matchInd$GenoBrNum],]
#######
pdDrop = pd[(!pd$BrNum %in% matchInd$MethBrNum) & (!pd$BrNum %in% missingGenotypes),] # 1 donors
missing_geno = snpCorMat[jaffelab::ss(rownames(snpCorMat),"_",1)%in% pdDrop$BrNum,]
#missing_geno_subset = missing_geno[ , colSums(missing_geno>0.7,na.rm=T)>0 | jaffelab::ss(colnames(missing_geno),"_") %in% jaffelab::ss(rownames(missing_geno),"_",1)  ]
#pdf("qc/pheatmap_missing_geno_subset_55SNP.pdf",h=15,w=15,onefile=FALSE)
#pheatmap(missing_geno_subset, 
#		cluster_rows=T, 
#		cluster_cols=T,
#		color=col.pal,
#		fontsize=16)
#dev.off()
#missing_geno_subset = cbind(GenoBrNum=rownames(missing_geno_subset), data.frame(missing_geno_subset,row.names=NULL) )


#table(pdDrop$predictedSex==pdDrop$Sex)
#######
snpCorWithin = cor(beta.genotypes.flipped, use="pairwise.complete.obs")
l = list()
listInd = 1
for (i in 1:ncol(snpCorWithin) ) {
	matchInd = which(snpCorWithin[i,]>0.7)
	if (length(matchInd) > 1) {
	subSnpCor = snpCorWithin[matchInd,matchInd]
	ids = ss(colnames(subSnpCor),"_")
		if (length(unique(ids))>1) { 
			l[[listInd]] = subSnpCor
			listInd = listInd+1
		}
	}
}
length(unique(l))
unique(l)
snpMismatchWithin =  snpCorWithin[rownames(snpCorWithin)%in% unique(c(sapply(unique(l),colnames))),colnames(snpCorWithin)%in%unique(c(sapply(unique(l),colnames)))]

pdf("qc/SupplementalFigure_pheatmap_snpMismatchWithin_methylation_55SNP.pdf",h=15,w=15,onefile=FALSE,useDingbats=FALSE)
pheatmap(snpMismatchWithin, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=30)
dev.off()


snpMismatchWithin = cbind(MethBrNum=rownames(snpMismatchWithin), data.frame(snpMismatchWithin,row.names=NULL) )

#table(pd[pd$BrNum%in%unique(colnames(snpMismatchWithin)),'predictedSex']==pd[pd$BrNum%in%unique(colnames(snpMismatchWithin)),'Sex'])
#2 sex swaps obsevered (hence 3 sample swaps previously left undetected!)
### Samples with matching BrNums, low correlation
l = list()
listInd = 1
br = as.character(unique(pd$BrNum))

for (i in 1:length(br)) {
	matchInd = which(pd$BrNum0 == br[i])
	if (length(matchInd) > 1) {
	subSnpCor = snpCorWithin[matchInd,matchInd]
	if (sum(subSnpCor < .7) > 0) {
		l[[listInd]] = subSnpCor
		listInd = listInd+1
		}
	}
}
length(l)
unique(l)
## None

### Annotation from master list
masterlist_summary = openxlsx::read.xlsx(xlsxFile='/dcl01/lieber/ajaffe/lab/brain_swap/brainSwap_master.xlsx')
mismatched_geno$Geno_Masterlist_Drop = mismatched_geno$GenoBrNum %in% masterlist_summary$Drop
mismatched_geno$Meth_Masterlist_Drop = mismatched_geno$MethBrNum %in% masterlist_summary$Drop
mismatched_geno$Either_Masterlist_Drop = mismatched_geno$Geno_Masterlist_Drop| mismatched_geno$Meth_Masterlist_Drop 

mismatched_geno$Geno_Masterlist_Flag = mismatched_geno$GenoBrNum %in% masterlist_summary$Flag
mismatched_geno$Geno_Masterlist_FlagTimes = masterlist_summary[match(mismatched_geno$GenoBrNum, masterlist_summary$Flag),'FlagTimes']
mismatched_geno$Meth_Masterlist_Flag = mismatched_geno$MethBrNum %in% masterlist_summary$Flag
mismatched_geno$Meth_Masterlist_FlagTimes = masterlist_summary[match(mismatched_geno$MethBrNum, masterlist_summary$Flag),'FlagTimes']
mismatched_geno$Either_MasterlistFlag = mismatched_geno$Geno_Masterlist_Flag| mismatched_geno$Meth_Masterlist_Flag 

mismatched_geno$AnyProblem = mismatched_geno$Either_Masterlist_Drop| mismatched_geno$Either_MasterlistFlag 

### Finalized decision
mismatched_geno$Decision = "Drop"
#mismatched_geno$Decision[mismatched_geno$MethBrNum %in% c("Br5083") ] <- "Keep" #DROP GENOTYPES mismatched

openxlsx::write.xlsx(list(mismatchedGeno=mismatched_geno, `missingGeno`=pdDrop, snpMismatchWithin=snpMismatchWithin),file='qc/n398_meth450k_genotype_problematic_samples.xlsx')