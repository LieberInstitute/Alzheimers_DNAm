### Age acceleration and cellular composition estimates

### Composition
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
library(ggplot2)
library(ggrepel)
library(jaffelab)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##
load('rdas/horvath_epigenetic_clock_output.rda')
pd2 = pd
rm(pd)
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda')

##### Cell Type over Development Plot #########
dat <- pd[,c('BrNum', 'ES', 'NPC','DA_NEURON','NeuN_pos','NeuN_neg','Region','Dx','DxOrdinal','Age','Sex','Race')]
dat = tidyr::gather(dat, key="CellType", value="Proportion", ES,NPC,DA_NEURON,NeuN_pos,NeuN_neg)

cell_type_dx_subset = ggplot(data=dat[dat$DxOrdinal!="Alz Drop",], aes(x=Region,y=Proportion,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`CellType`,scales='free',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Age Group', y = "Cell Type Proportion", title='Subset Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 

pdf('plots/cell_proportions_by_dx_by_region.pdf',height=10,width=20)
cell_type_dx_subset
dev.off()

pdf('plots/SupplementalFigure_estimated_proportion_of_neurons_by_dx_by_region.pdf',height=12,width=12,useDingbats=FALSE)
cell_type_dx_subset = ggplot(data=dat[dat$DxOrdinal!="Alz Drop" & dat$CellType=="NeuN_pos",], aes(x=Region,y=Proportion,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85),size=2) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Brain Region', y = "Estimated Proportion of Neurons") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
  theme(legend.position = c(.875,.95), 
		legend.background = element_rect(colour = "black"),
		legend.title=element_blank(),
		legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )
print(cell_type_dx_subset)
dev.off()

########## Getting statistics
library(tidyr)
library(broom)
library(lmerTest)
straighten <- function(..., fn = tidy){
  fits <- list(...)
  if (is.null(names(fits))) names(fits) <- character(length(fits))
  # If a fit isn't named, use the object name
  dots <- match.call(expand.dots = FALSE)$...
  obj_nms <- vapply(dots, deparse, character(1))
  names(fits)[names(fits) == ""] <- obj_nms[names(fits) == ""]
  
  purrr::map2(.x = fits,
              .y = names(fits),
              .f = function(x, n){
                data.frame(model = n, 
                           fn(x),
                           stringsAsFactors = FALSE)
              }) %>%
    dplyr::bind_rows()
}
dat$DxOrdinal <- as.numeric(dat$DxOrdinal)

###########---Subset model stats ----############
#### all samples
subset_NeuN_neg_int_fixed = lm(NeuN_neg ~ Dx + Region + Dx:Region+Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop",])
subset_NeuN_pos_int_fixed = lm(NeuN_pos ~ Dx + Region + Dx:Region+Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop",])
subset_NeuN_neg_noInt_fixed = lm(NeuN_neg ~ Dx + Region + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop",])
subset_NeuN_pos_noInt_fixed = lm(NeuN_pos ~ Dx + Region + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop",])

#### stratified
## DLPFC
subset_NeuN_neg_noInt_dlpfc = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='DLPFC',])
subset_NeuN_pos_noInt_dlpfc = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='DLPFC',])
summary(subset_NeuN_pos_noInt_dlpfc)
## ERC
subset_NeuN_neg_noInt_erc = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='ERC',])
subset_NeuN_pos_noInt_erc = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='ERC',])
summary(subset_NeuN_pos_noInt_erc)
## HIPPO
subset_NeuN_neg_noInt_hippo = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='HIPPO',])
subset_NeuN_pos_noInt_hippo = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='HIPPO',])
summary(subset_NeuN_pos_noInt_hippo)
## CRB
subset_NeuN_neg_noInt_crb = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='CRB',])
subset_NeuN_pos_noInt_crb = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='CRB',])
summary(subset_NeuN_pos_noInt_crb)

###
t.test(NeuN_pos~Dx, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='DLPFC',])
t.test(NeuN_pos~Dx, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='ERC',])
t.test(NeuN_pos~Dx, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='HIPPO',])
t.test(NeuN_pos~Dx, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='CRB',])


### Save model results
res = straighten(
		   subset_NeuN_neg_int_fixed,
		   subset_NeuN_pos_int_fixed, 
		   subset_NeuN_neg_noInt_fixed,
		   subset_NeuN_pos_noInt_fixed)
		   
res$Model = jaffelab::ss(res$model,"_",1)		 
res$NeuN = jaffelab::ss(res$model,"_",3)		   
res$InteractionTerm = ifelse(jaffelab::ss(res$model,"_",4)=="int",TRUE,FALSE)	

res$term = gsub("DxOrdinal","DxAlzheimer", res$term)	
res$Sig = ifelse(res$`p.value`<0.05,TRUE,FALSE)	
res$Direction = ifelse(sign(res$estimate)==1, "Up", "Down")
 
res = res[order(res$term,res$`p.value`), ]
write.csv(res, file='csvs/cell_composition_linearModel_results.csv',row.names=F)

##### Age acceleration #########

##
pd$DNAmAge <- datout$DNAmAge[match(pd$Sample_Name,pd2$Sample_Name)]
cor.test(pd$DNAmAge,pd$Age)
cor.test(pd$DNAmAge,pd$negControl_PC1)

## some tests stratifeid by region
pd_tmp = pd[pd$DxOrdinal!="Alz Drop", ]
cor.test(pd_tmp$DNAmAge[pd_tmp$Region=="DLPFC"],pd_tmp$Age[pd_tmp$Region=="DLPFC"])
cor.test(pd_tmp$DNAmAge[pd_tmp$Region=="ERC"],pd_tmp$Age[pd_tmp$Region=="ERC"])
cor.test(pd_tmp$DNAmAge[pd_tmp$Region=="HIPPO"],pd_tmp$Age[pd_tmp$Region=="HIPPO"])
cor.test(pd_tmp$DNAmAge[pd_tmp$Region=="CRB"],pd_tmp$Age[pd_tmp$Region=="CRB"]) 
##
pd_dlpfc = pd_tmp[pd_tmp$Region=="DLPFC",]
pd_dlpfc$AA = residuals(lm(DNAmAge ~ Age, data =pd_dlpfc ) )
t.test(AA~Dx, pd_dlpfc)
##
pd_erc = pd_tmp[pd_tmp$Region=="ERC",]
pd_erc$AA = residuals(lm(DNAmAge ~ Age, data =pd_erc ) )
t.test(AA~Dx, pd_erc)
##
pd_hippo = pd_tmp[pd_tmp$Region=="HIPPO",]
pd_hippo$AA = residuals(lm(DNAmAge ~ Age, data =pd_hippo ) )
t.test(AA~Dx, pd_hippo)
##
pd_crb = pd_tmp[pd_tmp$Region=="CRB",]
pd_crb$AA = residuals(lm(DNAmAge ~ Age, data =pd_crb ) )
t.test(AA~Dx, pd_crb)
##
pd_tmp2 = rbind(pd_dlpfc, pd_erc, pd_hippo, pd_crb)

chronAge_v_DNAmAge = ggplot(data=pd[pd$DxOrdinal!="Alz Drop", ], aes(x=Age,y=DNAmAge )) + 
geom_point(aes(col=`Dx`),size=2) +
geom_smooth(se=FALSE,method='lm') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Chronological Age', y = "DNAm Age") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom',legend.title=element_blank(),panel.spacing = unit(2, "lines")) + 
    scale_x_continuous(breaks=c(50, 70, 90))
ggsave(chronAge_v_DNAmAge,filename='plots/SupplementalFigure_DNAmAge_v_chronologicalAge_byRegion.pdf',height=9,width=12,useDingbats=FALSE)

subset_dx_AA = ggplot(data=pd_tmp2, aes(x=Dx,y=AA,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85),size=2) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Diagnosis', y = "Age Acceleration") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none',axis.title.x=element_blank() ) 
ggsave(subset_dx_AA,filename='plots/SupplementalFigure_dx_vs_age_acceleration_by_region_subset_stratified_analysis.pdf',height=10,width=10,useDingbats=FALSE)
