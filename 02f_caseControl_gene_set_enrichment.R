###### Test all-region models
library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')

load('rdas/caseControl_DMC_allRegion.rda')

library(org.Hs.eg.db)
library(limma)

flatten_annotation <- function (dat=allStats, annotation_col='within10kb_geneSymbol_gencode_hg38') {
missing <-  dat[,annotation_col]==""
ann.keep <- dat[!missing,c('Name',annotation_col)]
geneslist <- strsplit(ann.keep[,annotation_col], split = ";")
names(geneslist) <- ann.keep$Name
    flat <- data.frame(symbol = unlist(geneslist))
    flat$symbol <- as.character(flat$symbol)
    flat$cpg <- substr(rownames(flat), 1, 10)
    flat$alias <- limma::alias2SymbolTable(flat$symbol)
    eg <- toTable(org.Hs.egSYMBOL2EG)
    m <- match(flat$symbol, eg$symbol)
    flat$entrezid <- eg$gene_id[m]
    flat <- flat[!is.na(flat$entrezid), ]
    id <- paste(flat$cpg, flat$entrezid, sep = ".")
    d <- duplicated(id)
    flat.u <- flat[!d, ]
	return(flat.u)
	}
##
gomethGen= function(flat.u=flat.u, sig.cpg) {
    sig.cpg <- as.character(sig.cpg)
    sig.cpg <- sig.cpg[!is.na(sig.cpg)]
    all.cpg <- unique(flat.u$cpg)
    sig.cpg <- unique(sig.cpg)
    m1 <- match(flat.u$cpg, sig.cpg)
    eg.sig <- flat.u$entrezid[!is.na(m1)]
    eg.sig <- unique(eg.sig)
    m2 <- match(flat.u$cpg, all.cpg)
    eg.all <- flat.u$entrezid[!is.na(m2)]
    freq_genes <- table(eg.all)
    eg.universe <- names(freq_genes)
    test.de <- as.integer(eg.universe %in% eg.sig)
    sorted.eg.sig <- eg.universe[test.de == 1]
    out <- list(sig.eg = sorted.eg.sig, universe = eg.universe,
        freq = freq_genes, de = test.de)
pwf <- missMethyl:::.estimatePWF(D = test.de, bias = as.vector(freq_genes))
GOgst <- limma::goana(sorted.eg.sig, universe = eg.universe,
                prior.prob = pwf)
GOgst$FDR <- p.adjust(GOgst$P.DE, method = "BH")
GOgst = GOgst[order(GOgst$P.DE),]				
				
KEGGgst <- limma::kegga(sorted.eg.sig, universe = eg.universe,
                prior.prob = pwf)
KEGGgst$FDR <- p.adjust(KEGGgst$P.DE, method = "BH")
KEGGgst = KEGGgst[order(KEGGgst$P.DE),]				
goRes = list(GO=GOgst, KEGG=KEGGgst)
return(goRes)
}				

flat.u = flatten_annotation()
ResUnsplit = gomethGen(flat.u=flat.u,sig.cpg=allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05,'Name']) 
ResUnsplit = c(split(ResUnsplit[[1]],ResUnsplit[[1]]$Ont), ResUnsplit[2])
ResUnsplit = lapply(ResUnsplit, tibble::rownames_to_column,var="ID")
names(ResUnsplit)[1:3] = paste0("GO_",names(ResUnsplit)[1:3]) 

##
ResHyper = gomethGen(flat.u=flat.u,sig.cpg=allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05 & allStats$Primary_subset_mainEffect_logFC>0,'Name']) 
ResHyper = c(split(ResHyper[[1]],ResHyper[[1]]$Ont), ResHyper[2])
ResHyper = lapply(ResHyper, tibble::rownames_to_column,var="ID")
names(ResHyper)[1:3] = paste0("GO_",names(ResHyper)[1:3]) 

##
ResHypo = gomethGen(flat.u=flat.u,sig.cpg=allStats[allStats$Primary_subset_mainEffect_adj.P.Val<0.05 & allStats$Primary_subset_mainEffect_logFC<0,'Name']) 
ResHypo = c(split(ResHypo[[1]],ResHypo[[1]]$Ont), ResHypo[2])
ResHypo = lapply(ResHypo, tibble::rownames_to_column,var="ID")
names(ResHypo)[1:3] = paste0("GO_",names(ResHypo)[1:3]) 

##
ResInter = gomethGen(flat.u=flat.u,sig.cpg=allStats[allStats$Primary_subset_interactionEffect_adj.P.Val<0.05,'Name']) 
ResInter = c(split(ResInter[[1]],ResInter[[1]]$Ont), ResInter[2])
ResInter = lapply(ResInter, tibble::rownames_to_column,var="ID")
names(ResInter)[1:3] = paste0("GO_",names(ResInter)[1:3]) 


## Save results
openxlsx::write.xlsx(ResUnsplit, file='csvs/SupplementalTable_CrossRegion_GO_Results_Unsplit.xlsx')
openxlsx::write.xlsx(ResHyper, file='csvs/SupplementalTable_CrossRegion_GO_Results_Hypermethylated.xlsx')
openxlsx::write.xlsx(ResHypo, file='csvs/SupplementalTable_CrossRegion_GO_Results_Hypomethylated.xlsx')
openxlsx::write.xlsx(ResInter, file='csvs/SupplementalTable_RegionDependent_GO_Results.xlsx')

#pdf('plots/gene_differentialMethylation_CpG_bias.pdf')
#missMethyl:::.plotBias(D = test.de, bias = as.vector(freq_genes))
#dev.off()
   # prior.prob <- bias
   # o <- order(bias)
   # prior.prob[o] <- tricubeMovingAverage(D[o], span = 0.5)
   # prior.prob
