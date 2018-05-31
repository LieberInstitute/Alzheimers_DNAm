## Load libraries
library('limma')
library('edgeR')
library('jaffelab')
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')

#### Load methylation  ####
## Processed data
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')
pd_meth = pd

#### Load RNAseq  ####

## Find samples
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
load('rdas/caseControl_DMC_allRegion.rda')

main_sigProbe = allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05, 'Name']
main_goi  =lapply( strsplit(allStats[main_sigProbe,'within10kb_EnsId_gencode_hg38'], ";"), unique )

dropId = which(lengths(main_goi)==0)

foi = data.frame(Probe = main_sigProbe[-dropId], Gene = unlist(lapply(main_goi[-dropId],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "allRegion_Main"
foi$MethP = allStats[foi$Probe, 'Primary_subset_mainEffect_P.Value' ]
foi$MethLogFC = allStats[foi$Probe, 'Primary_subset_mainEffect_logFC' ]

foi_ss = foi[foi$Gene %in% results$gencodeID,]

meth_ss = bVals[as.character(foi_ss$Probe),pd_meth$Sample_Name]
geneRpkm_ss = geneRpkm[ unique(foi_ss$Gene), pd_rna$SAMPLE_ID]

dat = cbind(pd_rna, t(meth_ss), t(geneRpkm_ss) )
dat$Dx = factor(dat$Dx, levels=c("Control","Alzheimer") )
## Stats
full_res = list()
for (i in 1:nrow(foi_ss) ){

ensembl_gene = foi_ss$Gene[i]
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
res$GeneSymbol = geneMap[foi_ss$Gene[i],'Symbol']
res$Coefficient = gsub("mod","",rownames(res))
rownames(res) = NULL
res = res[res$Coefficient==cpg,]
full_res[[i]] = res 
cat('.')
}

full_res = do.call(rbind, full_res)
full_res$minCorrelation_pvalue = matrixStats::rowMins( as.matrix(full_res[,c("CRB_Pr(>|t|)", "ERC_Pr(>|t|)", "HIPPO_Pr(>|t|)", "DLPFC_Pr(>|t|)")]))				 
full_res$region_minP = apply( as.matrix(full_res[,c("CRB_Pr(>|t|)", "ERC_Pr(>|t|)", "HIPPO_Pr(>|t|)", "DLPFC_Pr(>|t|)")]),1, which.min)
full_res$region_minP = sapply(full_res$region_minP, function(x) c("CRB", "ERC", "HIPPO", "DLPFC")[x] )

save(full_res, file='rdas/DMP_beta_vs_gene_expression_mainEffect_Genes.rda')
write.csv(full_res, file='csvs/DMP_beta_vs_gene_expression_mainEffect_Genes.csv',row.names=F)

#### get number with possible association with DNAm ####
rownames(full_res) <- full_res$EnsemblGene
load('rdas/diffential_expression_statistics_for_DMP_Genes.rda',verbose=T)

nom_de = full_res[full_res$EnsemblGene %in% unique(foiMain$gencodeID[foiMain$minExprs_pvalue<0.05]),]
length(unique(nom_de[nom_de$minCorrelation_pvalue<0.05,'EnsemblGene']))

#### Now plotting these results ####
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))		

full_res_ss=nom_de[nom_de$minCorrelation_pvalue<0.05, ]
full_res_ss$min_Bonf = sapply(full_res_ss$minCorrelation_pvalue,p.adjust,method='bonf',n=length(unique(foiMain$Gene[foiMain$minExprs_pvalue<0.05])) )

write.csv(full_res_ss, file='csvs/SupplementalTable_top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.csv',row.names=F)			 

pdf('plots/SupplementalFigure_top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.pdf')
for (i in 1:nrow(full_res_ss)) {

ensembl_gene = full_res_ss$EnsemblGene[i]
cpg = as.character(full_res_ss$CpG[i])

pval = full_res_ss[i,'minCorrelation_pvalue']
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

a=ggplot(data=dat, aes_string(x=cpg,y=ensembl_gene) ) + geom_point(aes(col=Dx)) + facet_wrap(~Region,ncol=4) + labs(x=paste0( cpg," DNAm (%)"),
y=paste0( full_res_ss$GeneSymbol[i],"\n Adjusted RPKM"))+
geom_smooth(method='lm',se=F) +
scale_colour_manual(values=c("black","#ab1323" ) ) + 
theme(legend.position='none')
print(a) }
dev.off()

