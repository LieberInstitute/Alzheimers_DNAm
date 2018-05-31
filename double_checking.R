## Double checking stats changing
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
oldDat=bVals
##
load('/dcl01/lieber/ajaffe/Steve/Alz/Paper/rdas/cleanSamples_n377_processed_data_postfiltered.rda')
newDat=bVals

difVal = which(rowSums(oldDat[,colnames(newDat)]!=newDat)>0)

oldDat[difVal,colnames(newDat)] - newDat[difVal,]

difs = matrixStats::rowMaxs(oldDat[,colnames(newDat)] - newDat[,])
difs=difs[difs!=0]

######
pca = prcomp(t(bVals))
pca=pca$x
dat = cbind(pd,pca)
PC1_2 <- 

ggplot(data=dat[dat$DxOrdinal!="Alz Drop",] ,aes(x=PC3, y=Age, col=flagSample) ) +
	 geom_point(size=5)+ 
  theme(legend.position = c(.10,.15), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")
ggsave(PC1_2, file='qc/SupplementalFigure_alz_brainRegion_pca.pdf',height=8,width=12,useDingbats=FALSE )