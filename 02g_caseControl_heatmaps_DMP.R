## Alz heatmaps
library(pheatmap)
library(RColorBrewer)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## load data
setwd('/dcl01/lieber/ajaffe/Steve/Alz')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda',verbose=T)
##
pd$keepList[is.na(pd$keepList)] = TRUE
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") )


## load in stats
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')

allStats = allStats[order(allStats$`Primary_subset_interactionEffect_P.Value`),]
int_sigCpG = allStats[allStats$`Primary_subset_interactionEffect_adj.P.Val`<0.05,'Name']

## 1000 most significant interaction CpGs
int_sigCpG=int_sigCpG[1:(min(length(int_sigCpG),100))]

subsetIndex=which(pd$keepList)
mod <- model.matrix(~Dx+ negControl_PC1 + negControl_PC2 + Age+Sex + snpPC1 + Region + Dx:Region, data = pd[subsetIndex,])
mod = mod[,c(1:2,8:13,3:7)]
clean_dat = t( jaffelab::cleaningY(bVals[int_sigCpG, subsetIndex],mod=mod,P=8  ) ) 
			
			
##########		
pd_subset = pd[subsetIndex,]
annotation <- data.frame( Dx = pd_subset$Dx,
						  Region = pd_subset$Region,
#						  Race = pd_subset$Race,
#						  Sex = pd_subset$Sex,
						  stringsAsFactors=T)
rownames(annotation) <- pd_subset$Sample_Name
library(RColorBrewer)
ann_colors = list(Dx = c('red','blue') ,
				  Region = brewer.pal(8,"Set1")[3:6])
	#			  Race= brewer.pal(8,"Dark2")[1:2],
#				  Sex= c('black','grey') )
names(ann_colors[['Dx']])<- levels(annotation$Dx)
names(ann_colors[['Region']])<- levels(annotation$Region)
#names(ann_colors[['Sex']])<- levels(annotation$Sex)
#names(ann_colors[['Race']])<- levels(annotation$Race)

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/clustered_cleaned_bVals_interaction_CpGs_FDR05_top1000_subsetModel.pdf',height=20,width=30,onefile=F)
pheatmap::pheatmap(clean_dat,
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   treeheight_row=200,
				   treeheight_col=200,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = TRUE,
				   fontsize = 14,
				   scale='column')
dev.off()	 

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
ann_colors = list(Dx = c('y','red') )
#				  Region = brewer.pal(8,"Set1")[3:6],
#				  Race= brewer.pal(8,"Dark2")[1:2],
#				  Sex= c('black','grey') )
names(ann_colors[['Dx']])<- levels(annotation$Dx)
#names(ann_colors[['Region']])<- levels(annotation$Region)
#names(ann_colors[['Sex']])<- levels(annotation$Sex)
#names(ann_colors[['Race']])<- levels(annotation$Race)

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

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/clustered_cleaned_bVals_mainEffect_CpGs_FDR05_top100_subsetModel.pdf',height=20,width=20,onefile=F)
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

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/clustered_uncleaned_bVals_mainEffect_CpGs_FDR05_top1000_subsetModel.pdf',height=20,width=30,onefile=F)
pheatmap::pheatmap(t(bVals[main_sigCpG, subsetIndex]),
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   treeheight_row=200,
				   treeheight_col=200,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = FALSE,
				   fontsize = 14,
				   scale='none')
dev.off()

