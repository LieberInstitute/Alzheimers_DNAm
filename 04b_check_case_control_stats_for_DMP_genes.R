## Load libraries
library('limma')
library('edgeR')
library('jaffelab')
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
#### Load RNAseq  gene-level results ####
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_results_hg38.Rdata',verbose=T)

########### Methylation: allRegion stats###########
load('rdas/caseControl_DMC_allRegion.rda')

### Pull main dm probes as well as associated genes
main_sigProbe = allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05, 'Name']
main_goi  =lapply( strsplit(allStats[main_sigProbe,'within10kb_EnsId_gencode_hg38'], ";"), unique )

dropId = which(lengths(main_goi)==0)

foi = data.frame(Probe = main_sigProbe[-dropId], Gene = unlist(lapply(main_goi[-dropId],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "allRegion_Main"
foi$MethP = allStats[foi$Probe, 'Primary_subset_mainEffect_P.Value' ]
foi$MethLogFC = allStats[foi$Probe, 'Primary_subset_mainEffect_logFC' ]

## Pull rna-seq results
dlpfcResults = results[results$Region=="DLPFC",]
ercResults = results[results$Region=="ERC",]
cbResults = results[results$Region=="CB",]
hippoResults = results[results$Region=="HIPPO",]

foiMain = cbind(foi, 
				 dlpfcResults[match(foi$Gene, dlpfcResults$gencodeID ),c('gencodeID','Symbol', 'log2FC','pvalue','qvalue','meanExprs')],
				 ercResults[match(foi$Gene, ercResults$gencodeID ),c('log2FC','pvalue','qvalue','meanExprs')],
				 cbResults[match(foi$Gene, cbResults$gencodeID ),c('log2FC','pvalue','qvalue','meanExprs')],
				 hippoResults[match(foi$Gene, hippoResults$gencodeID ),c('log2FC','pvalue','qvalue','meanExprs')])
colnames(foiMain)[8:ncol(foiMain)]	  <- paste0( rep(c("DLPFC_","ERC_","CRB_","HIPPO_"),each=4), colnames(foiMain)[8:ncol(foiMain)]	) 
				 
foiMain=foiMain[!is.na(foiMain$gencodeID),]
foiMain$minExprs_pvalue = matrixStats::rowMins( as.matrix(foiMain[,c("DLPFC_pvalue", "ERC_pvalue", "CRB_pvalue", "HIPPO_pvalue")]))		
foiMain$region_minP = apply( as.matrix(foiMain[,c("DLPFC_pvalue", "ERC_pvalue", "CRB_pvalue", "HIPPO_pvalue")]),1, which.min)
foiMain$region_minP = sapply(foiMain$region_minP, function(x) c("DLPFC", "ERC", "CRB", "HIPPO")[x] )

foiMain$min_bonferroni_P =  sapply(foiMain$minExprs_pvalue,p.adjust,method='bonf',n=length(unique(foiMain$Gene)) )

table(foiMain$minExprs_pvalue[!duplicated(foiMain$Gene)]<0.05)
table(foiMain$min_bonferroni_P[!duplicated(foiMain$Gene)]<0.05)
unique(foiMain[foiMain$min_bonferroni_P<0.05,'Symbol'])

save(foiMain, file='rdas/diffential_expression_statistics_for_DMP_Genes.rda')
write.csv(foiMain,file='csvs/SupplementalTable_diffential_expression_statistics_for_DMP_Genes.csv',row.names=FALSE)
