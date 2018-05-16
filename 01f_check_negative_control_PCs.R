#qsub -l bluejay -l mf=10G,h_vmem=20G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH 02_check_negative_control_PCs.R
# Analysis of Negative Control PCs: Hippocampus
setwd('/dcl01/lieber/ajaffe/Steve/Alz/')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/RGset_n398.rda')

library(ggplot2)
library(ggrepel)
library(jaffelab)
library(minfi)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

#############
sampleNames(RGset) <- pd$Sample_Name
#Fix pd brain number
pd$BrNum[!is.na(pd$BrNum)]= paste0("Br", as.character(as.numeric(gsub("Br", "", pd$BrNum[!is.na(pd$BrNum)]) )) )
#Change Hippocampus to Hippo
pd$Sample_Plate<-pd$Sample_Group
#############
# PCA of negative control PCs
controlProbes = minfi:::.extractFromRGSet450k(RGset)
pca = prcomp(t(log2(rbind(controlProbes$greenControls$NEGATIVE,
	controlProbes$redControls$NEGATIVE)+1)))
negControlPCs =	pca$x[,1:4]
colnames(negControlPCs) = paste0("negControl_", colnames(negControlPCs))
pd = cbind(pd, negControlPCs)	
pd$negControlPC_VarExplained <- jaffelab::getPcaVars(pca)

#########				 
pd$Sentrix_ID <- factor(pd$Sentrix_ID, levels = unique(pd$Sentrix_ID[order(pd$Sample_Plate,pd$Sentrix_ID)]))
#########				 
a <- ggplot(data=pd, aes(x=Sample_Plate, y=negControl_PC1,fill=Sample_Plate) ) +
	 geom_boxplot() + 
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
	 labs(y=paste0("NegControl PC1"," (", as.character(signif(pd$negControlPC_VarExplained[1],3)),"%)" )) + 
	 theme(axis.text.x = element_text(angle = 45, size = 12, hjust = .95), axis.title.x=element_blank() )
b <- ggplot(data=pd, aes(x=Sentrix_ID, y=negControl_PC1,fill=Sample_Plate) ) +
	 geom_boxplot() + 
	 scale_colour_brewer(palette = "Set3") +
	 scale_fill_brewer(palette = "Set3") + 
	 labs(y="NegControl PC1") +
     theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size =8),
	 axis.title=element_blank())		
c <- ggplot(data=pd, aes(x=Sample_Plate, y=negControl_PC2,fill=Sample_Plate) ) +
	 geom_boxplot() + 
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
	 labs(y=paste0("NegControl PC2"," (", as.character(signif(pd$negControlPC_VarExplained[2],3)),"%)" )) +		theme(axis.text.x = element_text(angle = 45, size = 12, hjust = .95), axis.title.x=element_blank() )
d <- ggplot(data=pd, aes(x=Sentrix_ID, y=negControl_PC2,fill=Sample_Plate) ) +
	 geom_boxplot() + 
	 scale_colour_brewer(palette = "Set3") +
	 scale_fill_brewer(palette = "Set3") + 
	 labs(y=paste0("NegControl PC2"," (", as.character(signif(pd$negControlPC_VarExplained[2],3)),"%)" )) +
     theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size =8),
	 axis.title=element_blank())	
#########
library(gridExtra)
library(grid)
lay <- rbind(c(1,2,2),
             c(3,4,4) )			   
top4_negControl_PCs <- arrangeGrob(a, b, c,d, ncol=2,nrow=2, layout_matrix = lay) #generates g
ggsave(top4_negControl_PCs, file="qc/top4_negControl_PCs.pdf", height=8.5,width=11) 	 
#########
pd$Sentrix_Pos_RN <- gsub("C.*","",pd$Sentrix_Position)

e <- ggplot(data=pd, aes(x=Sentrix_Position, y=negControl_PC1, fill=Sample_Plate) ) +
	 geom_boxplot()  + facet_wrap(~Sample_Plate,scales='fixed',ncol=6)+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
		labs(y="NegControl PC1") + 
		theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size = 8), axis.title.x=element_blank() )
ggsave(e, file="qc/Plate_by_Position_PC1.pdf", height=8.5,width=11) 	 

f <- ggplot(data=pd, aes(x=Sentrix_Pos_RN, y=negControl_PC1, fill=Sample_Plate) ) +
	 geom_boxplot()  + facet_wrap(~Sample_Plate,scales='fixed',ncol=6)+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
		labs(y="NegControl PC1") + 
		theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size = 8), axis.title.x=element_blank() )
ggsave(f, file="qc/Plate_by_Position_RN_PC1.pdf", height=8.5,width=11)

g <- ggplot(data=pd, aes(x=Sentrix_Position, y=negControl_PC2, fill=Sample_Plate) ) +
	 geom_boxplot()  + facet_wrap(~Sample_Plate,scales='fixed',ncol=6)+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
		labs(y="NegControl PC2") + 
		theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size = 8), axis.title.x=element_blank() )
ggsave(g, file="qc/Plate_by_Position_PC2.pdf", height=8.5,width=11) 	 

h <- ggplot(data=pd, aes(x=Sentrix_Pos_RN, y=negControl_PC2, fill=Sample_Plate) ) +
	 geom_boxplot()  + facet_wrap(~Sample_Plate,scales='fixed',ncol=6)+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + 
		labs(y="NegControl PC2") + 
		theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=0.3, size = 8), axis.title.x=element_blank() )
ggsave(h, file="qc/Plate_by_Position_RN_PC2.pdf", height=8.5,width=11) ###

########################## PCs across first Four
PC1_2 <- ggplot(data=pd, aes(x=negControl_PC1, y=negControl_PC2, col=Sample_Plate) ) +
	 geom_point()+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +
				theme(legend.position="bottom")

PC2_3 <- ggplot(data=pd, aes(x=negControl_PC2, y=negControl_PC3, col=Sample_Plate) ) +
	 geom_point()+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +
				 geom_text_repel(data=dplyr::filter(pd, abs(negControl_PC3)>10), aes(label=Sample_Name),colour="black")		

PC3_4 <- ggplot(data=pd, aes(x=negControl_PC3, y=negControl_PC4, col=Sample_Plate) ) +
	 geom_point()+
	 scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +
				 geom_text_repel(data=dplyr::filter(pd, abs(negControl_PC4)>10), aes(label=Sample_Name),colour="black")		 

dat <- data.frame(VarExplained=pd[1:20,'negControlPC_VarExplained'],index = 1:20)
PC_Scree <- ggplot(data=dat, aes(x=index, y= VarExplained)) + geom_point() + labs(x="negControl PC Number", y = "Variance Explained (%)") + geom_vline(xintercept=4.5,linetype=3)
negControlPC_scatter <- arrangeGrob(PC1_2, PC2_3, PC3_4,PC_Scree, ncol=2,nrow=2) 

#generates g
ggsave(negControlPC_scatter, file="qc/negControlPC_scatter.pdf", height=8.5,width=11)				 
				 