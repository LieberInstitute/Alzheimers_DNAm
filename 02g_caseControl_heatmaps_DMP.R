## Alz heatmaps
library(pheatmap)
library(RColorBrewer)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## load data
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda',verbose=T)

## load in stats
load('rdas/caseControl_DMC_allRegion.rda')

####### Main effects

allStats = allStats[order(allStats$`Primary_subset_mainEffect_P.Value`),]
main_sigCpG = allStats[allStats$`Primary_subset_mainEffect_adj.P.Val`<0.05,'Name']
#main_sigCpG=main_sigCpG[1:(min(length(main_sigCpG),100))]

subsetIndex=which(pd$keepList)

mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region, data = pd[subsetIndex,])
clean_dat=	t( jaffelab::cleaningY(bVals[main_sigCpG, subsetIndex],mod=mod,P=2  ) ) 
			
##########		
pd_subset = pd[subsetIndex,]
annotation <- data.frame( Dx = pd_subset$Dx )
#						  Region = pd_subset$Region,
#						  Race = pd_subset$Race,
#						  Sex = pd_subset$Sex,
#						  stringsAsFactors=T)
rownames(annotation) <- pd_subset$Sample_Name
library(RColorBrewer)
ann_colors = list(Dx = c("black","#ab1323" ) )
#				  Region = brewer.pal(8,"Set1")[3:6],
#				  Race= brewer.pal(8,"Dark2")[1:2],
#				  Sex= c('black','grey') )
names(ann_colors[['Dx']])<- levels(annotation$Dx)

col_scaled_clean_dat = apply(clean_dat, 2, scale)
rownames(col_scaled_clean_dat) <- rownames(clean_dat)
## changing color scale for clarity
quantile(col_scaled_clean_dat,c(.01,.99)) #~1.5 both sidescol_scaled_clean_dat =
col_scaled_clean_dat[col_scaled_clean_dat< -2.5] = -2.5
col_scaled_clean_dat[col_scaled_clean_dat> 2.5] = 2.5

breaksList = seq(-2.5, 2.5, length.out =100)

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(as.dist( col_scaled_clean_dat ) )
mat_cluster_cols <- sort_hclust(mat_cluster_cols)

pdf('plots/Figure_clustered_cleaned_bVals_mainEffect_CpGs_FDR05_subsetModel.pdf',height=20,width=22,onefile=F,useDingbats=FALSE)
pheatmap::pheatmap(col_scaled_clean_dat,
				   breaks = breaksList,
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   treeheight_row=200,
				   treeheight_col=0,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = FALSE,
				   fontsize = 25,
				   scale='none')
dev.off()


####### Main effects region effect left in

allStats = allStats[order(allStats$`Primary_subset_mainEffect_P.Value`),]
main_sigCpG = allStats[allStats$`Primary_subset_mainEffect_adj.P.Val`<0.05,'Name']
#main_sigCpG=main_sigCpG[1:(min(length(main_sigCpG),100))]

subsetIndex=which(pd$keepList)

mod <- model.matrix(~Dx+ Region+negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 , data = pd[subsetIndex,])
clean_dat=	t( jaffelab::cleaningY(bVals[main_sigCpG, subsetIndex],mod=mod,P=5  ) ) 
			
##########		
pd_subset = pd[subsetIndex,]
annotation <- data.frame( Dx = pd_subset$Dx )
#						  Region = pd_subset$Region,
#						  Race = pd_subset$Race,
#						  Sex = pd_subset$Sex,
#						  stringsAsFactors=T)
rownames(annotation) <- pd_subset$Sample_Name
library(RColorBrewer)
ann_colors = list(Dx = c("black","#ab1323" ) )
#				  Region = brewer.pal(8,"Set1")[3:6],
#				  Race= brewer.pal(8,"Dark2")[1:2],
#				  Sex= c('black','grey') )
names(ann_colors[['Dx']])<- levels(annotation$Dx)

col_scaled_clean_dat = apply(clean_dat, 2, scale)
rownames(col_scaled_clean_dat) <- rownames(clean_dat)
## changing color scale for clarity
quantile(col_scaled_clean_dat,c(.01,.99)) #~1.5 both sidescol_scaled_clean_dat =
col_scaled_clean_dat[col_scaled_clean_dat< -2.5] = -2.5
col_scaled_clean_dat[col_scaled_clean_dat> 2.5] = 2.5

breaksList = seq(-2.5, 2.5, length.out =100)

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(as.dist( col_scaled_clean_dat ) )
mat_cluster_cols <- sort_hclust(mat_cluster_cols)

pdf('plots/regionEffectLeftIn_Figure_clustered_cleaned_bVals_mainEffect_CpGs_FDR05_subsetModel.pdf',height=20,width=22,onefile=F,useDingbats=FALSE)
pheatmap::pheatmap(col_scaled_clean_dat,
				   breaks = breaksList,
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   treeheight_row=200,
				   treeheight_col=0,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = FALSE,
				   fontsize = 25,
				   scale='none')
dev.off()
