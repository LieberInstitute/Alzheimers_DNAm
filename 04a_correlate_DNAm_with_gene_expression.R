## Load libraries
library('limma')
library('edgeR')
library('jaffelab')

#### Load methylation  ####
## Processed data
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
pd_meth = pd
pd_meth$Dx = factor(pd_meth$Dx, levels=c("Control","Alzheimer") )
pd_meth$Region = factor(pd_meth$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd_meth$keepList[is.na(pd_meth$keepList)] = TRUE
pd_meth$DxOrdinal = as.character(pd_meth$Dx)
pd_meth[!pd_meth$keepList,'DxOrdinal'] <- 'Alz Drop'
pd_meth[pd_meth$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'

#### Load RNAseq  ####

## Find samples
#pd_rna <- read.csv('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis/GSK_PhaseII_SampleInformation.csv', header = TRUE)
#
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_results_hg38.Rdata',verbose=T)
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_pd_hg38.Rdata')
pd_rna = pd 
pd_rna$Region = plyr::revalue(pd_rna$Region, c("CB"="CRB") )
pd_rna$Dx = plyr::revalue(pd_rna$Dx, c("NC"="Control","AD" = "Alzheimer") )
pd_rna$BrNum = paste0("Br", pd_rna$BRNum)

## Load gene data
#load('/dcl01/ajaffe/data/lab/libd_alzheimers/hg38/rpkmCounts_alzheimer_gsk_phaseII_n422.rda')
#save(pd, geneRpkm,geneMap, file='/dcl01/lieber/ajaffe/Steve/Alz/rnaseq_data/geneRpkm.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rnaseq_data/geneRpkm.rda')

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
geneRpkm = geneRpkm[results$gencodeID,pd_rna$SAMPLE_ID]

###### Run analysis #######
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')

main_sigProbe = allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05, 'Name']
main_goi  =lapply( strsplit(allStats[main_sigProbe,'UCSC_RefGene_Name'], ";"), unique )

dropId = which(lengths(main_goi)==0)

foi = data.frame(Probe = main_sigProbe[-dropId], Gene = unlist(lapply(main_goi[-dropId],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "allRegion_Main"
foi$MethP = allStats[foi$Probe, 'Primary_subset_mainEffect_P.Value' ]
foi$MethLogFC = allStats[foi$Probe, 'Primary_subset_mainEffect_logFC' ]

foi_ss = foi[foi$Gene %in% geneMap[results$gencodeID,]$Symbol,]

meth_ss = bVals[as.character(foi_ss$Probe),pd_meth$Sample_Name]
geneRpkm_ss = geneRpkm[ unique(geneMap[match(foi_ss$Gene, geneMap$Symbol),'gencodeID']), pd_rna$SAMPLE_ID]

dat = cbind(pd_rna, t(meth_ss), t(geneRpkm_ss) )
dat$Dx = factor(dat$Dx, levels=c("Control","Alzheimer") )
## Stats
full_res = list()
for (i in 1:nrow(foi_ss) ){

ensembl_gene = geneMap[match(foi_ss$Gene[i], geneMap$Symbol),'gencodeID']
cpg = as.character(foi_ss$Probe[i])

res = do.call("cbind", lapply(unique(dat$Region), function(reg) {
dat_region = dat[dat$Region==reg, ]

mod = cbind( dat_region[,c('age','Sex','Dx','Race','mitoRate','RIN','gene_Assigned_Percent')], dat_region[,cpg ] )
colnames(mod)[ncol(mod)] <- cpg
mod = model.matrix(~., mod)

res = summary(lm( log2(dat_region[,ensembl_gene]+1) ~ mod+0 ))
res=as.data.frame(res$coefficients)
colnames(res) <- paste0(reg,"_", colnames(res) )
res
}) )

res$Model = paste0(cpg,"_",ensembl_gene)
res$CpG = cpg
res$EnsemblGene= ensembl_gene
res$GeneSymbol = foi_ss$Gene[i]
res$Coefficient = gsub("mod","",rownames(res))
rownames(res) = NULL
res = res[res$Coefficient==cpg,]
full_res[[i]] = res 
cat('.')
}

full_res = do.call(rbind, full_res)
full_res$minCorrelation_pvalue = matrixStats::rowMins( as.matrix(full_res[,c("CRB_Pr(>|t|)", "ERC_Pr(>|t|)", "HIPPO_Pr(>|t|)", "DLPFC_Pr(>|t|)")]))				 
#full_res = full_res[order(full_res$"minCorrelation_pvalue"),]

save(full_res, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/DMP_beta_vs_gene_expression_mainEffect_Genes.rda')
write.csv(full_res, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/DMP_beta_vs_gene_expression_mainEffect_Genes.csv',row.names=F)

#### get number with possible association with DNAm ####
rownames(full_res) <- full_res$EnsemblGene
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/diffential_expression_statistics_for_DMP_Genes.rda',verbose=T)
nom_de = full_res[full_res$EnsemblGene %in% unique(foiMain$gencodeID[foiMain$minExprs_pvalue<0.05]),]
length(unique(nom_de[nom_de$minCorrelation_pvalue<0.05,'EnsemblGene']))

#### Now plotting these results ####
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))		

full_res=nom_de[nom_de$minCorrelation_pvalue<0.05, ]
write.csv(full_res, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.csv',row.names=F)			 
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.pdf')
for (i in 1:min(nrow(full_res),30) ){

ensembl_gene = full_res$EnsemblGene[i]
cpg = as.character(full_res$CpG[i])

pval = full_res[i,'minCorrelation_pvalue']
#est = full_res_cpg_effects[match(Model, full_res_cpg_effects$Model),'Estimate']

mod = cbind( dat[,cpg ], dat[,c('age','Sex','Dx','Race','mitoRate','RIN','gene_Assigned_Percent')]  )
colnames(mod)[1] <- cpg
mod = model.matrix(~., mod)

	cleanDat = lapply(unique(dat$Region), function(x) {
	keepRegion = which(dat$Region==x)
	dat_ss = dat[keepRegion, ]
	mod_ss=mod[keepRegion, ]
	cleanGeneRpkm = jaffelab::cleaningY(y= log2(geneRpkm_ss[ensembl_gene,dat_ss$SAMPLE_ID,drop=F]+1 ), mod=mod_ss,P=2)
	
	return(list(dat_ss,cleanGeneRpkm) )} )
	
pd_clean = do.call( "rbind", lapply(cleanDat, `[[`, 1) )
gene_clean = do.call( "cbind", lapply(cleanDat, `[[`, 2) )
dat = cbind(pd_clean, t(gene_clean) )

a=ggplot(data=dat, aes_string(x=cpg,y=ensembl_gene) ) + geom_point() + facet_wrap(~Region,ncol=4) + labs(x=paste0( cpg," DNAm (%)"),
y=paste0( full_res$GeneSymbol[i],"\n Adjusted RPKM"))+
#title = paste0(cpg, " vs. ", full_res$GeneSymbol[i], "\n", "min p=", as.character(signif(pval,3))) ) +
geom_smooth(method='lm',se=F)
print(a)
}
dev.off()

######################### stratified stuff
## Stats
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
### Pull DLPFC probes as well as associated genes
DLPFC_sigProbe = regionSpecific_mergedStats[regionSpecific_mergedStats$DLPFC_subset_NoAdj_adj.P.Val<0.10, 'Name']
DLPFC_goi  =lapply( strsplit(subsetRegionSig[DLPFC_sigProbe,'UCSC_RefGene_Name'], ";"), unique )

dropId = which(lengths(DLPFC_goi)==0)

foi = data.frame(Probe = DLPFC_sigProbe[-dropId], Gene = unlist(lapply(DLPFC_goi[-dropId],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "DLPFC"
foi$MethP = regionSpecific_mergedStats[foi$Probe, 'DLPFC_subset_NoAdj_P.Value' ]
foi$MethLogFC = regionSpecific_mergedStats[foi$Probe, 'DLPFC_subset_NoAdj_logFC' ]

dlpfcResults[match(foi$Gene, dlpfcResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs','gencodeID','NumTx')]


foi_ss = foi[foi$Gene %in% geneMap$Symbol,]

meth_ss = bVals[as.character(foi_ss$Probe),pd_meth$Sample_Name]
geneRpkm_ss = geneRpkm[match(foi_ss$Gene, geneMap$Symbol), pd_rna$SAMPLE_ID]

dat = cbind(pd_rna, t(meth_ss), t(geneRpkm_ss) )
dat$Dx = factor(dat$Dx, levels=c("Control","Alzheimer") )

full_res = list()
for (i in 1:nrow(foi_ss) ){

ensembl_gene = geneMap[match(foi_ss$Gene[i], geneMap$Symbol),'gencodeID']
cpg = as.character(foi_ss$Probe[i])


mod = cbind( dat[,c('age','Sex','Dx','Race','mitoRate','Region','RIN','gene_Assigned_Percent')], dat[,cpg ] )
colnames(mod)[ncol(mod)] <- cpg
mod = model.matrix(~., mod)

res = summary(lm( log2(dat[,ensembl_gene]+1) ~ mod+0 ))
res=as.data.frame(res$coefficients)
res$Model = paste0(cpg,"_",ensembl_gene)
res$CpG = cpg
res$EnsemblGene= ensembl_gene
res$GeneSymbol = foi_ss$Gene[i]
res$Coefficient = gsub("mod","",rownames(res))
rownames(res) = NULL
full_res[[i]] = res 
cat('.')
}
full_res = do.call(rbind, full_res)

full_res$CpG = jaffelab::ss(full_res$Model,"_",1)
full_res$Gene = jaffelab::ss(full_res$Model,"_",2)
full_res_cpg_effects = full_res[full_res$Coefficient==full_res$CpG,]
full_res_cpg_effects = full_res_cpg_effects[order(full_res_cpg_effects$"Pr(>|t|)"),]

save(full_res, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/DMP_beta_vs_gene_expression.rda')
write.csv(full_res_cpg_effects, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/DMP_beta_vs_gene_expression.csv',row.names=F)

### Now plotting these results ###
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/top_DMP_methylation_versus_corresponding_gene_expression.pdf')
for (i in 1:nrow(foi_ss) ){

ensembl_gene = geneMap[match(foi_ss$Gene[i], geneMap$Symbol),'gencodeID']
cpg = as.character(foi_ss$Probe[i])
Model=paste0(cpg,"_",ensembl_gene)
pval = full_res_cpg_effects[match(Model, full_res_cpg_effects$Model),'Pr(>|t|)']
est = full_res_cpg_effects[match(Model, full_res_cpg_effects$Model),'Estimate']

mod = cbind( dat[,cpg ], dat[,c('age','Sex','Dx','Race','mitoRate','Region','RIN','gene_Assigned_Percent')]  )
colnames(mod)[1] <- cpg
mod2 = model.matrix(~., mod)
cleanGene = jaffelab::cleaningY(y=t( log2(dat[,ensembl_gene,drop=FALSE]+1) ), mod=mod2,P=2)

plot_dat=cbind(mod,t(cleanGene) )
a=ggplot(data=plot_dat, aes_string(x=cpg,y=ensembl_gene) ) + geom_point() +labs(x=paste0( cpg,"\nDNA Methylation (%)"),
y=paste0( foi_ss$Gene[i],"\nAdjusted RPKM"),
title = paste0(cpg, " vs. ", foi_ss$Gene[i], "\n", "p=", as.character(signif(pval,3)), "\n","estimate=", as.character(signif(est,3)) )) +geom_smooth(method='lm',se=F)
print(a)
}
dev.off()
