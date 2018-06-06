library(minfi)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(jaffelab)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
## Region specific
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

### Pull in protein-coding genes for plotting
gencode_v25_GR38= rtracklayer::import(con = "/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
protein_coding_genes_gcV25 = unique( gencode_v25_GR38[gencode_v25_GR38$gene_type=="protein_coding",]$gene_name )

pd$Dx = dplyr::recode(pd$Dx, Control = "Control", Alzheimer = "Alzheimer's")
#### ERC Boxplots ####
regionSpecific_mergedStats = regionSpecific_mergedStats[order(regionSpecific_mergedStats$`ERC_subset_NoAdj_P.Value`),]
ERC_sigCpG = regionSpecific_mergedStats[regionSpecific_mergedStats$ERC_subset_NoAdj_adj.P.Val<0.1,'Name']

ercSubset_keep = which(pd$keepList & pd$Region=="ERC")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[ercSubset_keep,])
dat = cbind(pd[ercSubset_keep,],
			t( jaffelab::cleaningY(bVals[ERC_sigCpG,ercSubset_keep],mod=mod,P=2  ) ) )

### Create boxplots
pdf('plots/SupplementalFigure_best_ERC_subset_DMC_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (cpg_i in ERC_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, regionSpecific_mergedStats$Name)

pc_genes=unique(unlist(strsplit(regionSpecific_mergedStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )
custom_title = paste0( 
"ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
"\n ",fixName ) #custom title

a = ggplot(dat, aes_string(x = 'Dx', y = cpg_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(aes(col=`Dx`),width=0.3,size=2 ) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x="Diagnosis", title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black')) 
	
print(a)			
}
dev.off()

#### DLPFC Boxplots ####
regionSpecific_mergedStats = regionSpecific_mergedStats[order(regionSpecific_mergedStats$`DLPFC_subset_NoAdj_P.Value`),]
DLPFC_sigCpG = regionSpecific_mergedStats[regionSpecific_mergedStats$DLPFC_subset_NoAdj_adj.P.Val<0.1,'Name']

dlpfcSubset_keep = which(pd$keepList & pd$Region=="DLPFC")
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1, data = pd[dlpfcSubset_keep,])
dat = cbind(pd[dlpfcSubset_keep,],
			t( jaffelab::cleaningY(bVals[DLPFC_sigCpG, dlpfcSubset_keep],mod=mod,P=2  ) ) )

### Create boxplots
pdf('plots/SupplementalFigure_best_DLPFC_subset_DMC_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (cpg_i in DLPFC_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, regionSpecific_mergedStats$Name)

pc_genes=unique(unlist(strsplit(regionSpecific_mergedStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )
custom_title = paste0( 
"DLPFC p=",as.character(signif(regionSpecific_mergedStats[ii,'DLPFC_subset_NoAdj_P.Value'],3)), 
"\n ",fixName ) #custom title

a = ggplot(dat, aes_string(x = 'Dx', y = cpg_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(aes(col=`Dx`),width=0.3 ) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x="Diagnosis", title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none')  +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
	
print(a)			
}
dev.off()

#####################--------- All regions boxplots -------#############################
#### mainEffect Boxplots #####
load('rdas/caseControl_DMC_allRegion.rda')

allStats = allStats[order(allStats$`Primary_subset_mainEffect_P.Value`),]
main_sigCpG = allStats[allStats$`Primary_subset_mainEffect_adj.P.Val`<0.05,'Name']
main_sigCpG=main_sigCpG[1:(min(length(main_sigCpG),100))]

subsetIndex=which(pd$keepList)

mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex,])
dat = cbind(pd[subsetIndex,],
			t( jaffelab::cleaningY(bVals[main_sigCpG, subsetIndex],mod=mod,P=2  ) ) )
			
### Create boxplots
pdf('plots/Figure_best_mainEffect_subset_DMC_boxplots.pdf',height=10,width=12)
for (cpg_i in main_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, allStats$Name)
pc_genes=unique(unlist(strsplit(allStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )

custom_title = paste0( fixName ) #custom title


a = ggplot(dat, aes_string(x = 'Dx', y = cpg_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(aes(col=`Dx`),width=0.3 ) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x="Diagnosis", title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
		
print(a)			
}
dev.off()

######### Ordinal boxplots sensitivity example
allStats = allStats[order(allStats$`Primary_subset_mainEffect_P.Value`),]
main_sigCpG = allStats[allStats$`Primary_subset_mainEffect_adj.P.Val`<0.05,'Name']
main_sigCpG=main_sigCpG[1:(min(length(main_sigCpG),100))]

pd$DxOrdinal=factor(pd$DxOrdinal,levels=c('Control','Alz Drop','Alz Keep') )
mod <- model.matrix(~as.numeric(DxOrdinal)+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd)

dat = cbind(pd[],
			t( jaffelab::cleaningY(bVals[main_sigCpG, ],mod=mod,P=2  ) ) )

dat$DxOrdinal = plyr::revalue(dat$DxOrdinal, c('Alz Drop'='Asymptomatic\n Alzheimers', 'Alz Keep'='Symptomatic\n Alzheimers') )

			### Create boxplots
pdf('plots/Figure_best_mainEffect_subset_DMC_boxplots_ordinalModel.pdf',height=10,width=12)
for (cpg_i in main_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, allStats$Name)
pc_genes=unique(unlist(strsplit(allStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )

custom_title = paste0( fixName ) #custom title


a = ggplot(dat, aes_string(x = 'DxOrdinal', y = cpg_i)) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(width=0.3 ) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x="Diagnosis", title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30,colour='black',angle = 45, hjust = 1))
print(a)			
}
dev.off()

#### mainEffect Boxplots: regional effect left in #####
load('rdas/caseControl_DMC_allRegion.rda')

allStats = allStats[order(allStats$`Primary_subset_mainEffect_P.Value`),]
main_sigCpG = allStats[allStats$`Primary_subset_mainEffect_adj.P.Val`<0.05,'Name']
main_sigCpG=main_sigCpG[1:(min(length(main_sigCpG),100))]

subsetIndex=which(pd$keepList)

mod <- model.matrix(~Dx+ Region+negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 , data = pd[subsetIndex,])
dat = cbind(pd[subsetIndex,],
			t( jaffelab::cleaningY(bVals[main_sigCpG, subsetIndex],mod=mod,P=5  ) ) )
			
### Create boxplots
pdf('plots/best_mainEffect_subset_DMC_boxplots_regionEffectLeftIn.pdf',height=10,width=12)
for (cpg_i in main_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, allStats$Name)
pc_genes=unique(unlist(strsplit(allStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )

custom_title = paste0( fixName ) #custom title


a = ggplot(dat, aes_string(x = 'Region', y = cpg_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x='Diagnosis', title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
print(a)			
}
dev.off()

#### interaction Boxplots ####
allStats = allStats[order(allStats$`Primary_subset_interactionEffect_P.Value`),]
int_sigCpG = allStats[allStats$`Primary_subset_interactionEffect_adj.P.Val`<0.1,'Name']
int_sigCpG=int_sigCpG[1:(min(length(int_sigCpG),100))]

subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region + Dx:Region, data = pd[subsetIndex,])
mod = mod[,c(1:2,8:13,3:7)]
dat = cbind(pd[subsetIndex,],
			t( jaffelab::cleaningY(bVals[int_sigCpG, subsetIndex],mod=mod,P=8  ) ) )

### Create boxplots
pdf('plots/Figure_best_interactionEffect_subset_DMC_boxplots.pdf',height=10,width=12)
for (cpg_i in int_sigCpG) {

#Change column name for ggplot2 to work
ii=match(cpg_i, allStats$Name)
pc_genes=unique(unlist(strsplit(allStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]

fixName= paste(pc_genes, collapse = ", " )

custom_title = paste0( fixName ) #custom title

a = ggplot(dat, aes_string(x = 'Region', y = cpg_i, fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y=paste0(cpg_i, "\nDNAm level"),x='Diagnosis', title = custom_title) + 
#		scale_colour_brewer(palette = "Set1") +
#		scale_fill_brewer(palette = "Set1") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
print(a)			
}
dev.off()


######### OLD
##### DLPFC Boxplots ####
#mergedStats = mergedStats[order(mergedStats$`DLPFC_subset_P.Value`),]
#DLPFC_sigCpG = mergedStats[mergedStats$DLPFC_subset_adj.P.Val<0.1,'Name']
#
#### Create boxplots
#pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/best_DLPFC_subset_DMC_boxplots.pdf',height=11,width=8.5)
#for (cpg_i in DLPFC_sigCpG) {
#
#dat = cbind(pd,t(bVals[cpg_i,, drop=FALSE])  )
#
##Change column name for ggplot2 to work
#ii=match(cpg_i, mergedStats$Name)
#fixName= paste(unique(unlist(strsplit(mergedStats[ii,'UCSC_RefGene_Name'], ";"))), collapse = ", " )
#custom_title = paste0(cpg_i, 
#"\n DLPFC p=",as.character(signif(mergedStats[ii,'DLPFC_subset_P.Value'],3)), 
#"\n ",fixName ) #custom title
#
#a = ggplot(dat, aes_string(x = 'Region', y = cpg_i, fill='Dx')) +
#        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
#		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
#		labs(y=paste0(cpg_i, "% Methylation (beta)"), title = custom_title) + 
##		scale_colour_brewer(palette = "Set1") +
##		scale_fill_brewer(palette = "Set1") + 
#		scale_fill_manual(values=c("black","steelblue","#ab1323" ) ) + 
#		scale_colour_manual(values=c("black","steelblue","#ab1323" ) ) + 
#		theme(legend.position='bottom') 
#	
#print(a)			
#}
#dev.off()