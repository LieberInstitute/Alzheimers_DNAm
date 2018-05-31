## First determine potential background--it the intersection of all genes that could have been tagged by DNAm and all genes that could have been differentially expressed
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
##### GET METHYLATION GENES
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

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
######### STRINGdb
library(STRINGdb)
library(igraph)
geneList <- read.csv('csvs/SupplementalTable_top_DMP_methylation_versus_corresponding_gene_expression_nomDE_with_nomAssoc_mainEffect.csv')

###### Score 200
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=200, input_directory="" )
backgroundV = string_db$map( data.frame(gene=as.character(unique(GeneBackground))), "gene", removeUnmappedRows = TRUE, takeFirst=TRUE )
string_db$set_background(backgroundV$STRING_id)

V(background)$name

###### Score 400
geneList_proteinMapped <- string_db$map( data.frame(gene=as.character(unique(geneList$GeneSymbol[geneList$minCorrelation_pvalue<0.05]))), "gene", removeUnmappedRows = TRUE, takeFirst=TRUE )
string_db$get_ppi_enrichment(unique(geneList_proteinMapped$STRING_id))


degree(my.subgraph)==0
string_db$get_proteins()
string_db$score_threshold
degree(string_db$graph)==0

string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
backgroundV <- example1_mapped$STRING_id[1:2000]
## any score
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )

#string_db$set_background(geneList_proteinMapped$STRING_id ) # THIS IS CRITICAL FOR GETTING THE RIGHT PVALUE

geneList_proteinMapped <- string_db$map( data.frame(gene=unique(geneList$GeneSymbol)), "gene", removeUnmappedRows = TRUE )
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/string_ppi_networks_exprsDnaM_sig_score0.pdf')
string_db$plot_network( geneList_proteinMapped$STRING_id, add_link=FALSE, add_summary=TRUE)
## Fast greedy
clustersList = string_db$get_clusters(geneList_proteinMapped$STRING_id , algorithm="fastgreedy")
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]], add_link=FALSE, add_summary=TRUE)
}
dev.off()




string_db$get_ppi_enrichment(unique(geneList_proteinMapped$STRING_id))
string_db$get_ppi_enrichment(geneList_proteinMapped$STRING_id)


pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/string_ppi_networks_exprsDnaM_sig_score400.pdf')
string_db$plot_network( unique(geneList_proteinMapped$STRING_id), add_link=FALSE, add_summary=TRUE)
dev.off()

## Fast greedy
clustersList = string_db$get_clusters(geneList_proteinMapped$STRING_id , algorithm="fastgreedy")
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]], add_link=FALSE, add_summary=TRUE)
}

# get their neighbors
(V(graph)[neighbors(string_db$graph, geneList_proteinMapped$STRING_id[1])]$name)

neighbors(string_db$graph, geneList_proteinMapped$STRING_id[1])
neighbors(string_db$graph, geneList_proteinMapped$STRING_id[45])

L<-geneList_proteinMapped$STRING_id[geneList_proteinMapped$STRING_id%in% V(string_db$graph)$name]
a <- string_db$get_neighbors(L)

function (graph, v, mode = c("out", "in", "all", "total"))
{
    if (is.character(mode)) {
        mode <- igraph.match.arg(mode)
        mode <- switch(mode, out = 1, `in` = 2, all = 3, total = 3)
    }
    on.exit(.Call(C_R_igraph_finalizer))
    res <- .Call(C_R_igraph_neighbors, graph, as.igraph.vs(graph,
        v) - 1, as.numeric(mode))
    V(graph)[res + 1]
}

# get the subgraph of interest (your genes + their neighbors)
my.subgraph <- string_db$get_subnetwork(c(L, a))
my.subgraph <- string_db$get_subnetwork(c(L) )


trimmed_subgraph=delete.vertices(simplify(my.subgraph), degree(my.subgraph)==0)

## Build layout
l <- layout_with_fr(trimmed_subgraph)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf('plots/Figure_string_ppi_networks_subgraph_someInteraction_igraphed_score400.pdf',height=12,width=12)
plot(trimmed_subgraph, vertex.size=2, vertex.label=geneList_proteinMapped[match(V(trimmed_subgraph)$name, geneList_proteinMapped$STRING_id ),'gene'], vertex.label.dist=0.6, layout = l*2, vertex.label.font=3) 
dev.off()

png('plots/Figure_string_ppi_networks_subgraph_someInteraction_igraphed_score400.png',units="in",height=10,width=10,res=72)
plot(trimmed_subgraph, vertex.size=2, vertex.label=geneList_proteinMapped[match(V(trimmed_subgraph)$name, geneList_proteinMapped$STRING_id ),'gene'], vertex.label.dist=0.6, layout = l*2, vertex.label.font=3) 
dev.off()


# look how many genes and interactions you have
V(my.subgraph)
E(my.subgraph)
## 700 score
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=700, input_directory="" )
geneList_proteinMapped <- string_db$map( data.frame(gene=unique(geneList$GeneSymbol)), "gene", removeUnmappedRows = TRUE )
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/string_ppi_networks_exprsDnaM_sig_score700.pdf')
string_db$plot_network( geneList_proteinMapped$STRING_id, add_link=FALSE, add_summary=TRUE)
## Fast greedy
clustersList = string_db$get_clusters(geneList_proteinMapped$STRING_id , algorithm="fastgreedy")
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]], add_link=FALSE, add_summary=TRUE)
}
dev.off()


###########################
ppi_network=string_db$graph
 lambda = ppie.compLambda(igraph::degree(ppi_network, hitListSlice[hitListSlice>=0]), totalEdgeNum)


###
string_db$plot_ppi_enrichment( geneList_proteinMapped$STRING_id )


##
string_db$get_ppi_enrichment(geneList_proteinMapped$STRING_id)
string_db$get_ppi_enrichment(Alz.pathway.proteins)
Alz.pathway.proteins <- string_db$get_term_proteins('05010')$STRING_id

## KEGG Enrichment: would need to set background for this p-values to be right
enrichmentKEGG <- string_db$get_enrichment( geneList_proteinMapped$STRING_id, category = "KEGG", methodMT = "fdr", iea = TRUE )
enrichmentKEGG[enrichmentKEGG$pvalue_fdr<0.05,]

## GO Enrichment
enrichmentGO_bp <- string_db$get_enrichment( geneList_proteinMapped$STRING_id, category = "Process", methodMT = "fdr", iea = TRUE )