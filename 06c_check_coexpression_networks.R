### Check genes that replicated in Lunnon et al and are functionally implicated
load('rdas/DMP_beta_vs_gene_expression_mainEffect_Genes.rda')
lunnon<-read.csv('csvs/SupplementalTable_dmp_replicated_nomSig_dirConsistent_in_lunnon2014.csv')
lunnon[,c('ALL_main_Dx_P.Value_Lunnon2014')]

priority_res = read.csv('csvs/SupplementalTable_top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.csv')

subres = priority_res[priority_res$CpG%in%lunnon$Name,] 

subres = subres[order(subres$minCorrelation_pvalue),] 



##### GET METHYLATION GENES
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n377_processed_data_postfiltered.rda')
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

### drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/goldset_GencodeAnnotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
#7966
length(drop_hg38_unmappable) 

###
bVals <- bVals[-drop_hg38_unmappable, ] 					  
goldsetSub <- goldset[match(rownames(bVals),goldset$Name), ]					  
goldsetSub = plyr::rename(goldsetSub, c('predictedPos'='pos_hg38','pos'='pos_hg19','chr'='chr_hg19') )

DNAm_Genes = unique(unlist(strsplit(goldsetSub$within10kb_geneSymbol_gencode_hg38, split=';') ))
###### GET RNASEQ GENES
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_results_hg38.Rdata',verbose=T)
RNAseq_Genes = unique(results$Symbol)
RNAseq_Genes= RNAseq_Genes[RNAseq_Genes!=""]

GeneBackground = intersect(DNAm_Genes, RNAseq_Genes)
GeneForeground = as.character(unique(priority_res[priority_res$minCorrelation_pvalue<0.05,'GeneSymbol']) )
#library('org.Hs.eg.db')
#library('DOSE')
#library('ReactomePA')
#library(clusterProfiler)
#foreground = bitr(unique(subres[subres$minCorrelation_pvalue<0.05,'GeneSymbol']), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#background = bitr(GeneBackground, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#########
zhang_net_descriptions = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/Zhang_et_al_coexpression_networks.xlsx',sheet=2)

zhang_networks = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/Zhang_et_al_coexpression_networks.xlsx',sheet=4)
zhang_networks = zhang_networks[zhang_networks$Gene_Symbol %in% GeneBackground,]
zhang_networks = split(zhang_networks, zhang_networks$Module)
zhang_networks=zhang_networks[sapply(zhang_networks,nrow)>10]

tab = lapply(zhang_networks, function(x)  {

priority_and_network = length(intersect(GeneForeground,x$Gene_Symbol))
not_priority_and_network = length(x$Gene_Symbol) - priority_and_network
priority_and_not_network = length(GeneForeground) - priority_and_network
not_priority_and_not_network = sum(!(GeneBackground %in% x$Gene_Symbol) & !(GeneBackground %in% GeneForeground) )


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

#########
mostafavi_networks = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/Mostafavi_et_al_coexpression_networks.xlsx',sheet=1,startRow=4)
mostafavi_networks = mostafavi_networks[mostafavi_networks$Gene.Symbol %in% GeneBackground,]
mostafavi_networks = split(mostafavi_networks, mostafavi_networks$Module.ID)
mostafavi_networks=mostafavi_networks[sapply(mostafavi_networks,nrow)>10]

tab = lapply(mostafavi_networks, function(x)  {

priority_and_network = length(intersect(GeneForeground,x$Gene.Symbol))
not_priority_and_network = length(x$Gene.Symbol) - priority_and_network
priority_and_not_network = length(GeneForeground) - priority_and_network
not_priority_and_not_network = sum(!(GeneBackground %in% x$Gene.Symbol) & !(GeneBackground %in% GeneForeground) )


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
				
resSig$Bonf_Fisher = p.adjust(resSig$P_Fisher, method ='bonferroni', n = length(mostafavi_networks))
resSig$FDR_Fisher =  p.adjust(resSig$P_Fisher, method ='fdr', n = length(mostafavi_networks))
resSig = resSig[order(resSig$P_Fisher),]		


#ankPos=length(grep("ANK",unique(priority_res$GeneSymbol)))
#ankNeg=length(grep("ANK",unique(GeneBackground)))
#tab=matrix(c(ankPos,ankNeg,length(priority_res$GeneSymbol)-ankPos,length(GeneBackground)-ankNeg ),nrow=2 )
#fisher.test(tab)
