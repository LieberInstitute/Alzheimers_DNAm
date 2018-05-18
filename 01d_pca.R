#qsub -l bluejay,mf=30G,h_vmem=40G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 01d_pca.R

##PCA
library(ggplot2)
library(ggrepel)
library(jaffelab)
library(RColorBrewer)

theme_set(theme_bw(base_size=30) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##
load('rdas/processed_data_postfiltered.rda')

#### Modeling methylation changes
library(limma)

##
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(bVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
					  
### subset probes ###

# drop probes that do not map to hg38
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
drop_hg38_unmappable = which(!rownames(bVals) %in% goldset$Name)
length(drop_hg38_unmappable)
bVals <- bVals[-drop_hg38_unmappable, ] 

###
pca <- prcomp(t(bVals))
varExplained=jaffelab::getPcaVars(pca)
pca=pca$x
save(pd,pca,varExplained, file='rdas/pca_alz.rda')

#
colnames(pca) <- paste0("meth",colnames(pca) )


keep_list = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/keep_list.csv')
keep_list[!is.na(keep_list)] <- paste0("Br", as.character(as.numeric(keep_list[!is.na(keep_list)] )))

pd$keepList = ifelse(pd$Dx=="Control",NA,FALSE)
pd[pd$Region=="CRB" & pd$BrNum%in%keep_list$CB,'keepList' ] <- TRUE
pd[pd$Region=="DLPFC" & pd$BrNum%in%keep_list$DLPFC,'keepList' ] <- TRUE
pd[pd$Region=="ERC" & pd$BrNum%in%keep_list$ERC,'keepList' ] <- TRUE
pd[pd$Region=="HIPPO" & pd$BrNum%in%keep_list$Hippo,'keepList' ] <- TRUE


dat = cbind(pd,pca[,1:10] )
dat$Dx = factor(dat$Dx, levels=c("Control","Alzheimer") )
dat$Region = factor(dat$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
dat$keepList[is.na(dat$keepList)] = TRUE
dat$DxOrdinal = as.character(dat$Dx)
dat[!dat$keepList,'DxOrdinal'] <- 'Alz Drop'
dat[dat$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'

#> dat[dat$methPC1>0&dat$Region=="CRB",'BrNum']
#[1] "Br1909" "Br1615" "Br2257"

PC1_2 <- ggplot(data=dat[dat$DxOrdinal!="Alz Drop",] ,aes(x=methPC1, y=methPC2, col=Region) ) +
	 geom_point(size=5)+ 
	 labs(x=paste0("PC1 \n(",signif(varExplained[1],3) ,"% of variance)"),
		  y=paste0("PC2 \n(",signif(varExplained[2],3) ,"% of variance)")) +
  theme(legend.position = c(.10,.15), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")
ggsave(PC1_2, file='qc/SupplementalFigure_alz_brainRegion_pca.pdf',height=8,width=12,useDingbats=FALSE )

PC1_2 <- ggplot(data=dat[dat$DxOrdinal!="Alz Drop",] ,aes(x=methPC1, y=methPC2, col=Region) ) +
	 geom_boxplot(size=5)+ 
	 labs(x=paste0("PC1 \n(",signif(varExplained[1],3) ,"% of variance)"),
		  y=paste0("PC2 \n(",signif(varExplained[2],3) ,"% of variance)")) +
  theme(legend.position = c(.10,.15), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")
ggsave(PC1_2, file='qc/SupplementalFigure_alz_brainRegion_pca.pdf',height=8,width=12,useDingbats=FALSE )


PC2_3 <- ggplot(data=dat ,aes(x=methPC2, y=methPC3, col=Dx) ) +
	 geom_point(size=5)+ 
	 labs(x=paste0("PC1 \n(",signif(varExplained[1],3) ,"% of variance)"),
		  y=paste0("PC2 \n(",signif(varExplained[2],3) ,"% of variance)")) +
  theme(legend.position = c(.10,.15), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )  + scale_colour_brewer(palette="Dark2")

		
#
t.test(methPC3~Dx,data=dat)
ggsave(PC1_2, file='/dcl01/lieber/ajaffe/Steve/Alz/plots/alz_brainRegion_pca2_3.pdf',height=8,width=14 )
#ggsave(PC1_2, file='/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/brainseq_phase2/plots/development_versus_brainRegion_pca.png',height=8,width=14 )