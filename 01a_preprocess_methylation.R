#qsub -l mf=50G,h_vmem=90G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 01a_preprocess_methylation.R
library(minfi)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(jaffelab)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
###
setwd('/dcl01/lieber/ajaffe/Steve/Alz/')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/RGset_n398.rda')

#Fix pd brain number
pd$BrNum[!is.na(pd$BrNum)]= paste0("Br", as.character(as.numeric(gsub("Br", "", pd$BrNum[!is.na(pd$BrNum)]) )) )
###
sampleNames(RGset) = pd$Sample_Name

#############
# PCA of negative control PCs

controlProbes = minfi:::.extractFromRGSet450k(RGset)
pca = prcomp(t(log2(rbind(controlProbes$greenControls$NEGATIVE,
	controlProbes$redControls$NEGATIVE)+1)))
negControlPCs =	pca$x[,1:4]
colnames(negControlPCs) = paste0("negControl_", colnames(negControlPCs))
pd = cbind(pd, negControlPCs)	
pd$negControlPC_VarExplained <- jaffelab::getPcaVars(pca)

#############
## QC steps
Mset = mapToGenome(RGset)
qc = minfiQC(Mset)
qcMat = as.data.frame(qc$qc)
pd = cbind(pd, qcMat)

table(pd$predictedSex == pd$Sex)

############ QC Plots
qcPlot <- ggplot(data = pd, aes(x = mMed, y = uMed,col = factor(Region) ) )
qcPlot <- qcPlot + geom_point(shape=1)  +
			   labs(x = "Meth median intensity (log2)",
				    y = "Unmeth median intensity (log2)") + 	 
					scale_colour_brewer(palette = "Set1") +
					scale_fill_brewer(palette = "Set1") +
					theme(legend.position="bottom") + 
					geom_abline(intercept=10.5* 2,slope=-1,linetype=3) + 
					geom_text_repel(data=pd[which( (pd$mMed + pd$uMed)/2<9.5) ,], aes(label=Sample_Name),colour="black")
ggsave(qcPlot,filename='qc/qcPlot_Region.pdf')

### mean detection P-value plot
detP <- detectionP(RGset)
pd$mean_detection_p <- colMeans(detP)
pdf('qc/Mean_Detection_Pvalues.pdf')
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
hist(-log10(colMeans(detP)))
hist(-log10(colMeans(detP))[colMeans(detP)<0.1])
dev.off()
table(`detP_0.05`=colMeans(detP)>0.05,`detP_0.01`=colMeans(detP)>0.01 )

### qc plots
sn = paste0(pd$Region, ":", pd$Dx)
qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = sn, pdf = "qc/qcReport_RegionDx_n398.pdf")	

###
pd$Sentrix_Position = factor(pd$Sentrix_Position, 
	levels = unique(pd$Sentrix_Position))
oo = order(pd$Sentrix_ID, pd$Sentrix_Position)	
sn = paste0(pd$Region, ":", pd$Dx)

pdf(file='qc/qc_by_Sentrix_4regions2.pdf')
par(mar=c(7,5,2,2))
palette(brewer.pal(8,"Dark2"))
plot(qcMat$mMed[oo], bg = as.numeric(factor(pd$Sentrix_ID[oo])),
	pch = ifelse(pd$Sample_Group[oo] %in% paste("Plate", c(1,3,5)), 21, 22),
	ylab="Median log2 Meth Intensity", xaxt="n",xlab="",
	cex.axis=2,cex.lab=2,ylim=c(8,12.5))
abline(v=tapply(seq(along=pd$Sample_Group[oo]), pd$Sample_Group[oo], max)+0.5)
sIndexes=split0(pd$Sentrix_ID[oo])
axis(1,at=sapply(sIndexes[oo], mean), lab=names(sIndexes[oo]),las=2,cex.lab=1.4)
legend("bottomright", levels(factor(sn)),col=1:8, pch = 15,cex=1.3)	
plot(qcMat$uMed[oo], bg = as.numeric(factor(pd$Sentrix_ID[oo])),
	pch = ifelse(pd$Sample_Group[oo] %in% paste("Plate", c(1,3,5)), 21, 22),
	ylab="Median log2 Unmeth Intensity", xaxt="n",xlab="",
	cex.axis=2,cex.lab=2,ylim=c(8,12.5))
abline(v=tapply(seq(along=pd$Sample_Group[oo]), pd$Sample_Group[oo], max)+0.5)
sIndexes=split0(pd$Sentrix_ID[oo])
axis(1,at=sapply(sIndexes, mean), lab=names(sIndexes),las=2,cex.lab=1.4)
dev.off()

####### Predicted sex swaps
predictedSex_plot <- ggplot(data = pd[!is.na(pd$Sex),], aes(x = xMed, y = yMed) )
predictedSex_plot <- predictedSex_plot + geom_point(aes(colour = predictedSex) ) +
			   geom_point(data=pd[which(pd$predictedSex != pd$Sex),],
               pch=21, fill=NA, size=4, stroke=1, colour = 'black')	 +
			   labs(x = "X chr, median total intensity (log2)",
				    y = "Y chr, median total intensity (log2)",
				    title = "Predicted Sex by Methylation Intensity") +
				 geom_text_repel(data=pd[which(pd$predictedSex != pd$Sex),], aes(label=Sample_Name),colour="black")+
	 scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +
				theme(legend.position="bottom")
ggsave(predictedSex_plot,filename='qc/PredictedSex_by_MethylationIntensity.pdf') 

#### flag samples
pd$flagSample = cut(pd$mMed +pd$uMed, c(10,21,23,30))
levels(pd$flagSample) = c("LowQual", "CheckQual", "HighQual")


pdf("qc/medianIntensityPlot_gsk_n398.pdf")
library(RColorBrewer)
palette(brewer.pal(8,"Paired"))
plot(pd$mMed, pd$uMed, pch=21, bg=as.numeric(factor(sn)),
	xlab="Median Meth Intensity (log2)",
	ylab="Median Unmeth Intensity (log2)",
	main = "Intesity Plot", cex.axis=1.8, cex.lab=1.4, cex.main=1.5)
legend("topleft", levels(factor(sn)),col=1:8, pch = 15,cex=1.3)	
points(pd$uMed~ pd$mMed, subset = pd$flagSample=="LowQual", pch=0, cex=2)
points(pd$uMed~ pd$mMed, subset = pd$flagSample=="CheckQual", pch=2, cex=2)
legend("bottomright", c("LowQual", "CheckQual"), pch =c(0,2))
dev.off()

########## Dropping poor samples
table(pd$flagSample)
drop_index = which(pd$flagSample == "LowQual" )
pd = pd[-drop_index, ]
RGset = RGset[ ,-drop_index]
Mset = Mset[ ,-drop_index ]

############### Normalization
 Mset_SQN <- preprocessQuantile(RGset, 
	fixOutliers = TRUE, #Quantile normalization
	removeBadSamples = FALSE, 
	quantileNormalize = TRUE, 
	stratified = TRUE, 
	mergeManifest = TRUE, 
	sex = NULL)
	
### Get counts for predicting cell type	
load('/users/ajaffe/Lieber/Projects/450k/ECD2014/devMeth450k/rdas/cellComp_estimates_cellLines_NeuNs.rda')
counts = minfi:::projectCellType(getBeta(Mset_SQN[rownames(coefs),]), coefs)	
pd = cbind(pd, as.data.frame(counts[match(pd$Sample_Name, rownames(counts)),]) )	

save(Mset_SQN, pd, file = 'rdas/MsetSQN_prefiltered_n394.rda')
save(Mset, pd, file = 'rdas/MsetRaw_prefiltered_n394.rda')
save(RGset, pd, file = 'rdas/RGset_prefiltered_n394.rda')

################ visualise what the data looks like before and after normalisation
pdf('qc/beta_distribution_pre_post_SQN.pdf')
par(mfrow=c(1,2))
#density plot for raw data
densityPlot(RGset, main="Raw", legend=FALSE)
#density plot for normalized data
densityPlot(getBeta(Mset_SQN), main="SQN Normalized", legend=FALSE)
dev.off()

####################### Filtering
###
# remove any probes that have failed in one or more samples
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detectionP(RGset)
detP <- detP[match(featureNames(Mset_SQN),rownames(detP)),] 
failed <- detP > 0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
keep <- (rowMeans(failed)<0.05)  #removes probes that fail in over 5 percent of samples
table(keep) #5% rule "An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform"
Mset_SQN <- Mset_SQN[keep,] # dropping failed probes 

###
# Drop probes containing common SNP at SBE and target CpG, 1% MAF filter
Mset_SQN <- addSnpInfo(Mset_SQN)
orig=nrow(Mset_SQN)
orig

Mset_SQN <- dropLociWithSnps(Mset_SQN, snps=c("SBE","CpG"), maf=0.01)
new=nrow(Mset_SQN)
new

orig-new
#
####
##Drop probes with crossreactive behavior. 
##published in Chen et al. 2013
##http://www.tandfonline.com/doi/full/10.4161/epi.23470?scroll=top&needAccess=true
##Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray
xReactiveProbes <- read.csv(file="/users/ssemick/Tox_Expr/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
drop_xreactive <- !(featureNames(Mset_SQN) %in% xReactiveProbes$TargetID)
table(drop_xreactive)
Mset_SQN <- Mset_SQN[drop_xreactive,] 

### 
# Drop probes that map to the sex chromosomes
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
drop_sexChr <- !(featureNames(Mset_SQN) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(drop_sexChr)
Mset_SQN <- Mset_SQN[drop_sexChr,]

###
# filter people

keepIndex = pd$predictedSex==pd$Sex #keep people that do not have sex swaps
table(keepIndex)
Mset_SQN = Mset_SQN[,which(keepIndex)]
pd = pd[which(keepIndex), ]

###
# extract data on remaining people for diff methylation analysis
mVals <- getM(Mset_SQN)
bVals <- getBeta(Mset_SQN)

###
pdf('qc/Distribution_M_vals_and_Beta_vals_PostSQN.pdf')
par(mfrow=c(1,2))
densityPlot(bVals, main="Beta values", 
            legend=FALSE, xlab="Beta values")
densityPlot(mVals, main="M-values", 
            legend=FALSE, xlab="M-values")
dev.off()

######################
#More EDA via PCA Exploration
oo = order(matrixStats::rowSds(bVals),decreasing=TRUE)[1:100000]
pca = prcomp(t(bVals[oo,]))
save(pca, file = 'rdas/pca_100k_no_probe_filter.rda')
pcaVars = getPcaVars(pca)

PCs = pca$x[,1:15]
pd = cbind(as.data.frame(pd), PCs, pcaVars ) 

#################################
#Saving things out for meQTL analysis
save(Mset_SQN, pd, file = "rdas/MSet_SQN_postfiltered.rda")
save(pd, mVals, bVals, file = "rdas/processed_data_postfiltered.rda")