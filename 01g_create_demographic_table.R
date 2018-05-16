library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda',verbose=T)
##
pd$keepList[is.na(pd$keepList)] = TRUE
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") )


######## Get APOE4 risk allele based on RT qPCR
#apoe_qPCR <- read.csv('/dcl01/lieber/ajaffe/Steve/Alz/APOE_qPCR_genotypes_from_RT_031318.csv')
#pd[,c('rs7412','rs429358')] <- NA
#pd[,c('rs7412','rs429358')] <- apoe_qPCR[match(pd$BrNum, apoe_qPCR$`Br.`),c("rs7412", "rs429358")]
#table(rs429358=pd$rs429358,rs7412=pd$rs7412,useNA='ifany')
#
###
#apoe_risk_map <- read.csv('/dcl01/lieber/ajaffe/Steve/Alz/apoe_genotypes_from_SNP.csv')
#apoe_risk_map$rs7412_rs429358 = paste0(apoe_risk_map$rs7412,';',apoe_risk_map$rs429358)
#pd$APOE = apoe_risk_map[match( paste0(pd$rs7412,";",pd$rs429358), apoe_risk_map$rs7412_rs429358),'APOE']
#pd$APOE4_Dosage = sapply(strsplit( as.character(pd$APOE),"\\/" ), function(x) {sum(x=="E4")})
#save(pd, bVals,mVals,file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
fisher.test(table(pd$APOE4_Dosage[!duplicated(pd$BrNum)&pd$keepList], pd$Dx[!duplicated(pd$BrNum)&pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="CRB" & pd$keepList], pd$Dx[pd$Region=="CRB" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="DLPFC" & pd$keepList], pd$Dx[pd$Region=="DPLFC" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="HIPPO" & pd$keepList], pd$Dx[pd$Region=="HIPPO" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="ERC" & pd$keepList], pd$Dx[pd$Region=="ERC" & pd$keepList]))

#table(rowSums(table(pd$BrNum,pd$DxOrdinal)==0)!=2)

######### Pull based on SNP chip
#APOE epsilon4 risk allele status not possible to determine without imputation of rs429358 (not currently in our data)

### Load genotype data for samples passing QC
#load('/dcl01/lieber/ajaffe/Steve/Alz/genotype_data/Alz_meth450k_Genotypes_n380.rda')
#
#APOE SNPs: rs7412 rs429358
#### Extract SNPs of interest
#APOE_SNPS= paste0(snpMap$CHR,":",snpMap$POS) %in% c('19:45411941', '19:45412079') 
#ALZ_SNPS = snpMap[snpMap$name%in% c('rs429358','rs7412', 'rs2075650','rs11136000','rs3851179'),]
##rs2075650 COUNTED: "G"
##rs3851179 COUNTED: "T"
##rs7412 COUNTED: "T"
#
#ALZ_GENOTYPES = t(snp[rownames(ALZ_SNPS),])
#ALZ_GENOTYPES[match(pd$BrNum,rownames(ALZ_GENOTYPES) ),]
#colnames(ALZ_GENOTYPES) <- gsub( ":.*","", colnames(ALZ_GENOTYPES))
#
####
#pd = cbind(pd, ALZ_GENOTYPES[match(pd$BrNum,rownames(ALZ_GENOTYPES) ),])
#save(pd, bVals,mVals,file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
#######


##
library(LIBDpheno)
library(tableone)
pd = cbind(pd, toxicant[[1]][match(brnum(pd$BrNum),toxicant[[1]]$brnumerical),])
pd$APOE4_Dosage = factor(pd$APOE4_Dosage)
t1 <- CreateTableOne(vars = c("Age","Race","Sex","Sample_Group",'manner_of_death','cerad','braak','APOE4_Dosage'), strata = c("DxOrdinal","Region"), data = pd)
table1 <- print(t1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
write.csv(table1, file = "/dcl01/lieber/ajaffe/Steve/Alz/csvs/summary_of_sample_characteristics_table1.csv")


######
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
######
dat = pd[!duplicated(pd$BrNum),]
a <- ggplot(data=dat, aes(x=braak, y=Age) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2)		 
b <- ggplot(data=dat, aes(x=cerad, y=Age) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2) + theme(axis.text.x = element_text(angle = 45, hjust = .95, vjust=1, size = 12))
c <- ggplot(data=dat, aes(x=DxOrdinal, y=Age) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2)+ theme(axis.text.x = element_text(angle = 45, hjust = .95, vjust=1, size = 12))
d <- ggplot(data=dat, aes(x=DxOrdinal, y=brain_weight_gram) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2)		+ theme(axis.text.x = element_text(angle = 45, hjust = .95, vjust=1, size = 12))
e <- ggplot(data=dat, aes(x=DxOrdinal, y=as.numeric(braak)) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2)		+ theme(axis.text.x = element_text(angle = 45, hjust = .95, vjust=1, size = 12))
f <- ggplot(data=dat, aes(x=DxOrdinal, y=as.numeric(cerad)) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.2)		+ theme(axis.text.x = element_text(angle = 45, hjust = .95, vjust=1, size = 12))
library(gridExtra)
library(grid)
demo_boxplots <- arrangeGrob(a, b, c,d,e,f, ncol=3,nrow=2) 
ggsave(demo_boxplots, file="/dcl01/lieber/ajaffe/Steve/Alz/plots/demographic_boxplots.pdf", height=8,width=14) 	 	

### Get list of brnums to get CERAD +Braak data for
pd$braak_or_cerad_missing = ifelse(is.na(pd$braak) | is.na(pd$cerad), "Yes", "No" )
write.csv( pd[!duplicated(pd$BrNum),c('BrNum','braak_or_cerad_missing','source','braak','cerad','DxOrdinal')], file='/dcl01/lieber/ajaffe/Steve/Alz/Alz_98donors_neuropath_data_available_to_SS.csv',row.names=F )

## Random t-tests of interest
t.test(Age~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) 
t.test(brain_weight_gram~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) 
t.test(bmi_calculated~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) #
			