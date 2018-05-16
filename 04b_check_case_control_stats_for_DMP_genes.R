## Load libraries
library('limma')
library('edgeR')
library('jaffelab')

#### Load RNAseq  gene-level results ####
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_results_hg38.Rdata',verbose=T)

########### Methylation: Regional stats###########
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')

### Pull DLPFC probes as well as associated genes
DLPFC_sigProbe = regionSpecific_mergedStats[regionSpecific_mergedStats$DLPFC_subset_NoAdj_adj.P.Val<0.05, 'Name']
DLPFC_goi  =lapply( strsplit(regionSpecific_mergedStats[DLPFC_sigProbe,'UCSC_RefGene_Name'], ";"), unique )

#dropId = which(lengths(DLPFC_goi)==0)

foi = data.frame(Probe = DLPFC_sigProbe, Gene = unlist(lapply(DLPFC_goi,paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "DLPFC"
foi$MethP = regionSpecific_mergedStats[foi$Probe, 'DLPFC_subset_NoAdj_P.Value' ]
foi$MethLogFC = regionSpecific_mergedStats[foi$Probe, 'DLPFC_subset_NoAdj_logFC' ]
## Pull rna-seq results
dlpfcResults = results[results$Region=="ERC",]
foiDlpfc = cbind(foi, dlpfcResults[match(foi$Gene, dlpfcResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs','gencodeID','NumTx')])
foiDlpfc = foiDlpfc[order(foiDlpfc$pvalue),]

### Pull ERC probes as well as associated genes
ERC_sigProbe = regionSpecific_mergedStats[regionSpecific_mergedStats$ERC_subset_NoAdj_adj.P.Val<0.05, 'Name']
ERC_goi  =lapply( strsplit(regionSpecific_mergedStats[ERC_sigProbe,'UCSC_RefGene_Name'], ";"), unique )


foi = data.frame(Probe = ERC_sigProbe[], Gene = unlist(lapply(ERC_goi[],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "ERC"
foi$MethP = regionSpecific_mergedStats[foi$Probe, 'ERC_subset_NoAdj_P.Value' ]
foi$MethLogFC = regionSpecific_mergedStats[foi$Probe, 'ERC_subset_NoAdj_logFC' ]
##
ercResults = results[results$Region=="ERC",]
foiErc = cbind(foi, ercResults[match(foi$Gene, ercResults$Symbol) ,c('log2FC','pvalue','qvalue','meanExprs','gencodeID','NumTx')])
foiErc = foiErc[order(foiErc$pvalue),]
foiErc = foiErc[!is.na(foiErc$gencodeID),]

###
regional_foi = rbind(foiDlpfc,foiErc)
regional_foi = regional_foi[order(regional_foi$pvalue),]

########### Methylation: allRegion stats###########
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')

### Pull main dm probes as well as associated genes
main_sigProbe = allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05, 'Name']
main_goi  =lapply( strsplit(allStats[main_sigProbe,'UCSC_RefGene_Name'], ";"), unique )

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
				 dlpfcResults[match(foi$Gene, dlpfcResults$Symbol ),c('gencodeID', 'log2FC','pvalue','qvalue','meanExprs')],
				 ercResults[match(foi$Gene, ercResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')],
				 cbResults[match(foi$Gene, cbResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')],
				 hippoResults[match(foi$Gene, hippoResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')])
colnames(foiMain)[7:ncol(foiMain)]	  <- paste0( rep(c("DLPFC_","ERC_","CRB_","HIPPO_"),each=4), colnames(foiMain)[7:ncol(foiMain)]	) 
				 
foiMain=foiMain[!is.na(foiMain$gencodeID),]
foiMain$minExprs_pvalue = matrixStats::rowMins( as.matrix(foiMain[,c("DLPFC_pvalue", "ERC_pvalue", "CRB_pvalue", "HIPPO_pvalue")]))		
foiMain$minExprs_qvalue = matrixStats::rowMins( as.matrix(foiMain[,c("DLPFC_qvalue", "ERC_qvalue", "CRB_qvalue", "HIPPO_qvalue")]))				 
		 
foiMain = foiMain[order(foiMain$minExprs_pvalue),]

table(foiMain$minExprs_pvalue[!duplicated(foiMain$Gene)]<0.05)
table(foiMain$minExprs_qvalue[!duplicated(foiMain$Gene)]<0.05)

library(tidyr)
library(dplyr)
dat = foiMain[,c('Probe','Gene','MethP','MethLogFC','DLPFC_pvalue','ERC_pvalue','CRB_pvalue','HIPPO_pvalue','DLPFC_log2FC','ERC_log2FC', 'CRB_log2FC','HIPPO_log2FC')]
dat 

####
int_sigProbe = allStats[allStats$Primary_subset_interactionEffect_adj.P.Val<0.05, 'Name']
int_goi  =lapply( strsplit(allStats[int_sigProbe,'UCSC_RefGene_Name'], ";"), unique )

dropId = which(lengths(int_goi)==0)

foi = data.frame(Probe = int_sigProbe[-dropId], Gene = unlist(lapply(int_goi[-dropId],paste,collapse=';') ),stringsAsFactors=FALSE )

s <- strsplit(foi$Gene, split = ";")

foi = data.frame(Probe = rep(foi$Probe, lengths(s) ), Gene = unlist(s),stringsAsFactors=FALSE )
foi$sigRegion = "allRegion_Interaction"
foi$MethP = allStats[foi$Probe, 'Primary_subset_interactionEffect_P.Value' ]
###
foiInt = cbind(foi, 
				 dlpfcResults[match(foi$Gene, dlpfcResults$Symbol ),c('gencodeID', 'log2FC','pvalue','qvalue','meanExprs')],
				 ercResults[match(foi$Gene, ercResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')],
				 cbResults[match(foi$Gene, cbResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')],
				 hippoResults[match(foi$Gene, hippoResults$Symbol ),c('log2FC','pvalue','qvalue','meanExprs')])
colnames(foiInt)[6:ncol(foiInt)]	  <- paste0( rep(c("DLPFC_","ERC_","CRB_","HIPPO_"),each=4), colnames(foiInt)[6:ncol(foiInt)]	) 

foiInt=foiInt[!is.na(foiInt$gencodeID),]
foiInt$minExprs_pvalue = matrixStats::rowMins( as.matrix(foiInt[,c("DLPFC_pvalue", "ERC_pvalue", "CRB_pvalue", "HIPPO_pvalue")]))				 
foiInt = foiInt[order(foiInt$minExprs_pvalue),]
foiInt$MethLogFC = NA
#foiInt = foiInt[,colnames(foiMain)]

table(foiInt$minExprs_pvalue[!duplicated(foiInt$Gene)]<0.05)
### Merge all results then save
#allRegionFoi = rbind(foiMain, foiInt)
#col_reorder = c(grep("pvalue", colnames(allRegionFoi) ,value=T),
#				grep("log2FC", colnames(allRegionFoi) ,value=T),
#				grep("meanExprs", colnames(allRegionFoi) ,value=T) )
#col_reorder = c(colnames(allRegionFoi)[!colnames(allRegionFoi) %in% col_reorder], col_reorder)				
#allRegionFoi = allRegionFoi[,col_reorder]

save(regional_foi, foiMain, foiInt, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/diffential_expression_statistics_for_DMP_Genes.rda')

openxlsx::write.xlsx(list(regionSpecific=regional_foi, foiMain=foiMain, foiInt = foiInt), file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/diffential_expression_statistics_for_DMP_Genes.xlsx')
##########
table(foiMain$minExprs_pvalue[!duplicated(foiMain$Gene)]<0.05)
table(foiMain$minExprs_qvalue[!duplicated(foiMain$Gene)]<0.05)

library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

dat = foiMain[,c('Probe','Gene','MethP','MethLogFC','DLPFC_pvalue','ERC_pvalue','CRB_pvalue','HIPPO_pvalue','DLPFC_log2FC','ERC_log2FC', 'CRB_log2FC','HIPPO_log2FC')]
dat=gather(dat,key='Region_Stat',value='value',5:12 ) %>% separate("Region_Stat",into=c('region','statistic'),sep="_" )
dat = spread(dat,statistic,value)
dat$NominalSig = ifelse(dat$pvalue<0.05,"P<0.05","P>0.05")
dat$region=factor(dat$region,levels=c('DLPFC',"ERC","HIPPO",'CRB') )
ala = ggplot(data=dat, aes(x=-log10(MethP)*sign(MethLogFC),y=-log10(pvalue)*sign(log2FC) ) ) + facet_wrap(~region,ncol=4) + geom_point(aes(col=NominalSig))+ scale_colour_manual(values = c("black", "grey")) + geom_hline(yintercept=-1*-log10(0.05), linetype='dashed', color='steelblue', size=1)+ geom_hline(yintercept=-log10(0.05), linetype='dashed', color='steelblue', size=1) + labs(x='Signed -log10(MethP)', y="Signed -log10(ExprsP)")
ggsave(ala, filename = '/dcl01/lieber/ajaffe/Steve/Alz/plots/differential_methylation_vs_differential_expression.pdf',height=8,width=8 )
####3