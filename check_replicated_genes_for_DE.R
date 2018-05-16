### Check genes that replicated in Lunnon et al and are functionally implicated
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/DMP_beta_vs_gene_expression_mainEffect_Genes.rda')
lunnon<-read.csv('/dcl01/lieber/ajaffe/Steve/Alz/csvs/dmp_replicated_nomSig_dirConsistent_in_lunnon2014.csv')
lunnon[,c('ALL_main_Dx_P.Value_Lunnon2014')]

priority_res = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/csvs/top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.csv')

subres = priority_res[priority_res$CpG%in%lunnon$Name,] 

subres = subres[order(subres$minCorrelation_pvalue),] 



##### GET METHYLATION GENES
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  
goldsetSub = plyr::rename(goldsetSub, c('predictedPos'='pos_hg38','pos'='pos_hg19','chr'='chr_hg19') )
goldsetSub = goldsetSub[,c('chr_hg19','pos_hg19','chr_hg38','pos_hg38', intersect(colnames(goldsetSub), colnames(ann450kSub)) )]		

DNAm_Genes = unique(unlist(strsplit(goldsetSub$UCSC_RefGene_Name, split=';') ))
###### GET RNASEQ GENES
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_results_hg38.Rdata',verbose=T)
RNAseq_Genes = unique(results$Symbol)
RNAseq_Genes= RNAseq_Genes[RNAseq_Genes!=""]

GeneBackground = intersect(DNAm_Genes, RNAseq_Genes)

library('org.Hs.eg.db')
library('DOSE')
library('ReactomePA')

foreground = bitr(unique(priority_res[priority_res$minCorrelation_pvalue<0.05,'GeneSymbol']), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
background = bitr(GeneBackground, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


ggo <- enrichGO(gene     = foreground$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
			   keyType  = 'ENTREZID',
               readable = TRUE,
			   universe=background$ENTREZID)
			   
compareBP <- compareCluster(geneCluster = list(test=foreground$ENTREZID), universe=background$ENTREZID, ont='BP', fun = "enrichGO", OrgDb    = org.Hs.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.20, readable = TRUE)

compareMF <- compareCluster(geneCluster = list(test=foreground$ENTREZID), universe=background$ENTREZID, ont='MF', fun = "enrichGO", OrgDb    = org.Hs.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.20, readable = TRUE)

compareCC <- compareCluster(geneCluster = list(test=foreground$ENTREZID), universe=background$ENTREZID, ont='CC', fun = "enrichGO", OrgDb    = org.Hs.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.2, readable = TRUE)

compareReactome <- compareCluster(geneCluster = list(test=foreground$ENTREZID), universe=background$	ENTREZID, fun = "enrichPathway", qvalueCutoff = 0.2, pvalueCutoff = 0.2, readable = TRUE)

compareKegg = compareCluster(list(test=foreground$ENTREZID), fun = "enrichKEGG",organism = "human", 
                             universe = background$ENTREZID, qvalueCutoff = 0.2, pvalueCutoff = 0.20)

compareDO   = compareCluster(list(test=foreground$ENTREZID), fun = "enrichDO", universe = background$ENTREZID,
							 qvalueCutoff = 0.2, pvalueCutoff = 0.20, readable = TRUE)
			   			   
library("pathview")
dat = subres[subres$minCorrelation_pvalue<0.05,'ERC_Estimate'][!duplicated(subres$GeneSymbol)]
names(dat) = foreground$ENTREZID

hsa05010 <- pathview(gene.data  = dat,
                     pathway.id = "hsa05010",
                     species    = "hsa")

#########
zhang_net_descriptions = zhang_networks = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/Zhang_et_al_coexpression_networks.xlsx',sheet=2)

zhang_networks = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/Zhang_et_al_coexpression_networks.xlsx',sheet=4)
zhang_networks = zhang_networks[zhang_networks$Gene_Symbol %in% background$SYMBOL,]
zhang_networks = split(zhang_networks, zhang_networks$Module)
zhang_networks=zhang_networks[sapply(zhang_networks,nrow)>10]

tab = lapply(zhang_networks, function(x)  {

priority_and_network = length(intersect(foreground$SYMBOL,x$Gene_Symbol))
not_priority_and_network = length(x$Gene_Symbol) - priority_and_network
priority_and_not_network = length(foreground$SYMBOL) - priority_and_network
not_priority_and_not_network = sum(!(background$SYMBOL %in% x$Gene_Symbol) & !(background$SYMBOL %in% foreground$SYMBOL) )


tab = matrix(c(priority_and_network, #pos,pos
			 not_priority_and_network, #neg,pos
			 priority_and_not_network, #pos,neg
			 not_priority_and_not_network), #neg,neg
	  nr = 2, byrow=TRUE)

})

res = lapply(tab, function(tab) {try(fisher.test(tab) )} )

resSig=res[sapply(res,class)=="htest"]
resSig = data.frame(Module=names(resSig), 
				OR_Fisher = sapply(resSig,function(x) x$estimate ),
				P_Fisher=sapply(resSig,function(x) x$p),
				nGenesEnrich = sapply(tab, function(x) x[1,1]),
				nGenesNetwork = sapply(tab, function(x) x[1,1]+x[1,2]  ) )
				
resSig$Bonf_Fisher = p.adjust(resSig$P_Fisher, method ='bonferroni', n = length(zhang_networks))
resSig$FDR_Fisher =  p.adjust(resSig$P_Fisher, method ='fdr', n = length(zhang_networks))
resSig = cbind(resSig, zhang_net_descriptions[match(resSig$Module, zhang_net_descriptions$Module),])	

resSig = resSig[order(resSig$P_Fisher),]			
write.csv(resSig,'/dcl01/lieber/ajaffe/Steve/Alz/csvs/179prioritizedGene_enrichment_with_Zhang_Coexpression_networks.csv', row.names=FALSE)