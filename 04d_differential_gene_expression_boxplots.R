library(ggplot2)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

### Differential gene expression boxplots
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/diffential_expression_statistics_for_DMP_Genes.rda')
foiMain = foiMain[order(foiMain$MethP),]

### load gene expression data
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_pd_hg38.Rdata')
pd_rna = pd 
pd_rna$Region = plyr::revalue(pd_rna$Region, c("CB"="CRB") )
pd_rna$Dx = plyr::revalue(pd_rna$Dx, c("NC"="Control","AD" = "Alzheimer") )
pd_rna$BrNum = paste0("Br", pd_rna$BRNum)

## Load gene data
#load('/dcl01/ajaffe/data/lab/libd_alzheimers/hg38/rpkmCounts_alzheimer_gsk_phaseII_n422.rda')
#save(pd, geneRpkm,geneMap, file='/dcl01/lieber/ajaffe/Steve/Alz/rnaseq_data/geneRpkm.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rnaseq_data/geneRpkm.rda')


cleanDat = lapply(unique(pd_rna$Region), function(x) {
pd_ss = pd_rna[pd_rna$Region==x, ]

mod=model.matrix(~ Dx + RIN + age + Sex + Race + mitoRate + gene_Assigned_Percent,data=pd_ss)
cleanGeneRpkm = jaffelab::cleaningY(y= log2(geneRpkm[foiMain$gencodeID,pd_ss$SAMPLE_ID]+1 ), mod=mod,P=2)

return(list(pd_ss,cleanGeneRpkm) )})

pd_clean = do.call( "rbind", lapply(cleanDat, `[[`, 1) )
gene_clean = do.call( "cbind", lapply(cleanDat, `[[`, 2) )
dat = cbind(pd_clean, t(gene_clean) )

###### MAKE BOXPLOTS
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/differential_expression_boxplots_for_mainEffect_DMP_genes.pdf',height=10,width=12)
for (gene_i in unique(foiMain$gencodeID[foiMain$minExprs_pvalue<1e-2]) ) {

#Change column name for ggplot2 to work

#custom_title = paste0(foiMain[match(gene_i, foiMain$gencodeID),'Gene'], "\n",gene_i ) #custom title
custom_title = paste0( foiMain[match(gene_i, foiMain$gencodeID),'Gene'] ) #custom title

a = ggplot(dat, aes_string(x = 'Region', y = gene_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y=paste0( "Adjusted RPKM"),x='Diagnosis', title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
#foiMain[match(gene_i, foiMain$gencodeID),'Gene'],
print(a)			
}
dev.off()