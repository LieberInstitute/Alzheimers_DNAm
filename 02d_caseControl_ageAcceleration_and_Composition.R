### Age acceleration and cellular composition estimates

### Composition
## script for brainseq phase2 methylation results
library(ggplot2)
library(ggrepel)
library(jaffelab)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/horvath_epigenetic_clock_output.rda')
pd2 = pd
rm(pd)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/cleanSamples_n380_processed_data_postfiltered.rda')
#pd$ageGroup = cut(pd$Age, breaks = c(-0.5,0,0.6,10,20,50,100))
#levels(pd$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")
#pd$ageStage = ifelse(pd$Age<0, "Fetal", ifelse(pd$Age>17, "Adult", NA ))
pd$keepList[is.na(pd$keepList)] = TRUE
pd$Dx = factor(pd$Dx, levels=c("Control","Alzheimer") )
pd$Region = factor(pd$Region, levels=c("CRB","DLPFC","HIPPO","ERC") )
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") )

##### Cell Type over Development Plot #########
dat <- pd[,c('BrNum', 'ES', 'NPC','DA_NEURON','NeuN_pos','NeuN_neg','Region','Dx','DxOrdinal','Age','Sex','Race')]
dat = tidyr::gather(dat, key="CellType", value="Proportion", ES,NPC,DA_NEURON,NeuN_pos,NeuN_neg)

cell_type_dx_full = ggplot(data=dat, aes(x=Region,y=Proportion,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`CellType`,scales='free',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Age Group', y = "Cell Type Proportion", title='Full Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 

cell_type_dx_subset = ggplot(data=dat[dat$DxOrdinal!="Alz Drop",], aes(x=Region,y=Proportion,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`CellType`,scales='free',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Age Group', y = "Cell Type Proportion", title='Subset Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 

cell_type_dx_ordinal = ggplot(data=dat, aes(x=Region,y=Proportion,fill=DxOrdinal )) + 
geom_point(aes(col=`DxOrdinal`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`CellType`,scales='free',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Age Group', y = "Cell Type Proportion", title='Ordinal Model') + 
		scale_fill_manual(values=c("black","steelblue","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","steelblue","#ab1323" ) ) + 
		theme(legend.position='bottom') 

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/cell_proportions_by_dx_by_region.pdf',height=10,width=20)
cell_type_dx_full
cell_type_dx_subset
cell_type_dx_ordinal
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/estimated_proportion_of_neurons_by_dx_by_region.pdf',height=8,width=12)
cell_type_dx_subset = ggplot(data=dat[dat$DxOrdinal!="Alz Drop" & dat$CellType=="NeuN_pos",], aes(x=Region,y=Proportion,fill=Dx )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Brain Region', y = "Estimated Proportion of Neurons") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 
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
## ERC
subset_NeuN_neg_noInt_erc = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='ERC',])
subset_NeuN_pos_noInt_erc = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='ERC',])
## HIPPO
subset_NeuN_neg_noInt_hippo = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='HIPPO',])
subset_NeuN_pos_noInt_hippo = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='HIPPO',])
## CRB
subset_NeuN_neg_noInt_crb = lm(NeuN_neg ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='CRB',])
subset_NeuN_pos_noInt_crb = lm(NeuN_pos ~ Dx + Age + Sex + snpPC1, data = pd[pd$DxOrdinal!="Alz Drop" & pd$Region=='CRB',])

##
summary(subset_NeuN_pos_noInt_dlpfc)
summary(subset_NeuN_pos_noInt_erc)
summary(subset_NeuN_pos_noInt_hippo)
summary(subset_NeuN_pos_noInt_crb)


###########---Full model stats ----############
#### fixed effects
full_NeuN_neg_int_fixed = lm(Proportion ~ Dx + Region + Dx:Region +Age , data = dat[dat$CellType=="NeuN_neg",])
full_NeuN_pos_int_fixed = lm(Proportion ~ Dx + Region + Dx:Region +Age , data = dat[dat$CellType=="NeuN_pos",])
full_NeuN_neg_noInt_fixed = lm(Proportion ~ Dx + Region+Age , data = dat[dat$CellType=="NeuN_neg",])
full_NeuN_pos_noInt_fixed = lm(Proportion ~ Dx + Region+Age , data = dat[dat$CellType=="NeuN_pos",])
## mixed effects 
full_NeuN_neg_int_mixed = lmer(Proportion ~ Dx + Region + Dx:Region+ (1 | BrNum), data = dat[dat$CellType=="NeuN_neg",], REML = T) 
full_NeuN_pos_int_mixed = lmer(Proportion ~ Dx + Region + Dx:Region+ (1 | BrNum), data = dat[dat$CellType=="NeuN_pos",], REML = T)  
full_NeuN_pos_noInt_mixed = lmer(Proportion ~ Dx + Region + (1 | BrNum), data = dat[dat$CellType=="NeuN_neg",], REML = T)  
full_NeuN_neg_noInt_mixed = lmer(Proportion ~ Dx + Region + (1 | BrNum), data = dat[dat$CellType=="NeuN_pos",], REML = T)  


###########---Ordinal model stats ----############
#### stats
ordinal_NeuN_neg_int_fixed = lm(Proportion ~ DxOrdinal + Region + DxOrdinal:Region+Age, data = dat[dat$CellType=="NeuN_neg",])
ordinal_NeuN_pos_int_fixed = lm(Proportion ~ DxOrdinal + Region + DxOrdinal:Region+Age, data = dat[dat$CellType=="NeuN_pos",])
ordinal_NeuN_neg_noInt_fixed =lm(Proportion ~ DxOrdinal + Region+Age, data = dat[dat$CellType=="NeuN_neg",])
ordinal_NeuN_pos_noInt_fixed =lm(Proportion ~ DxOrdinal + Region+Age, data = dat[dat$CellType=="NeuN_pos",])
## mixed effects 
ordinal_NeuN_neg_int_mixed = lmer(Proportion ~ DxOrdinal + Region + DxOrdinal:Region+ (1 | BrNum), data = dat[dat$CellType=="NeuN_neg",], REML = T)  
ordinal_NeuN_pos_int_mixed = lmer(Proportion ~ DxOrdinal + Region + DxOrdinal:Region+ (1 | BrNum), data = dat[dat$CellType=="NeuN_pos",], REML = T)  
ordinal_NeuN_neg_noInt_mixed = lmer(Proportion ~ DxOrdinal + Region + (1 | BrNum), data = dat[dat$CellType=="NeuN_neg",], REML = T)  
ordinal_NeuN_pos_noInt_mixed = lmer(Proportion ~ DxOrdinal + Region + (1 | BrNum), data = dat[dat$CellType=="NeuN_pos",], REML = T)  

### Save model results
res = straighten(full_NeuN_neg_int_fixed,
		   full_NeuN_pos_int_fixed,
		   full_NeuN_neg_noInt_fixed,
		   full_NeuN_neg_noInt_fixed,
		   subset_NeuN_neg_int_fixed,
		   subset_NeuN_pos_int_fixed, 
		   subset_NeuN_neg_noInt_fixed,
		   subset_NeuN_pos_noInt_fixed,
		   ordinal_NeuN_neg_int_fixed, 
		   ordinal_NeuN_pos_int_fixed,
           ordinal_NeuN_neg_noInt_fixed, 
           ordinal_NeuN_pos_noInt_fixed)
res$Model = jaffelab::ss(res$model,"_",1)		 
res$NeuN = jaffelab::ss(res$model,"_",3)		   
res$InteractionTerm = ifelse(jaffelab::ss(res$model,"_",4)=="int",TRUE,FALSE)	

res$term = gsub("DxOrdinal","DxAlzheimer", res$term)	
res$Sig = ifelse(res$`p.value`<0.05,TRUE,FALSE)	
res$Direction = ifelse(sign(res$estimate)==1, "Up", "Down")
 
res = res[order(res$term,res$`p.value`), ]
write.csv(res, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/cell_composition_linearModel_results.csv',row.names=F)

##### Age acceleration #########
pd$DxOrdinal = as.character(pd$Dx)
pd[!pd$keepList,'DxOrdinal'] <- 'Alz Drop'
pd[pd$DxOrdinal=="Alzheimer",'DxOrdinal'] <- 'Alz Keep'
pd$DxOrdinal= factor(pd$DxOrdinal, levels=c("Control", "Alz Drop", "Alz Keep") )

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
geom_point(aes(col=`Dx`),) +
geom_smooth(se=FALSE,method='lm') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Chronological Age', y = "DNAm Age") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 
ggsave(chronAge_v_DNAmAge,filename='/dcl01/lieber/ajaffe/Steve/Alz/plots/DNAmAge_v_chronologicalAge_byRegion.pdf',height=10,width=10)

subset_dx_AA = ggplot(data=pd_tmp2, aes(x=Dx,y=AA )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age Acceleration") + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') 
ggsave(subset_dx_AA,filename='/dcl01/lieber/ajaffe/Steve/Alz/plots/dx_vs_age_acceleration_by_region_subset_stratified_analysis.pdf',height=10,width=10)


##
summary(lm(DNAmAge ~ Age + Region + Age:Region, pd))
pd$AA_Linear <- residuals(lm(DNAmAge ~ Age + Region + Age:Region, pd))
##
full_dx_AA = ggplot(data=pd, aes(x=Dx,y=AA_Linear )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age Acceleration", title='Full Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 
		
subset_dx_AA = ggplot(data=pd[pd$DxOrdinal!="Alz Drop", ], aes(x=Dx,y=AA_Linear )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age Acceleration", title='Subset Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 

ordinal_dx_AA = ggplot(data=pd, aes(x=DxOrdinal,y=AA_Linear )) + 
geom_point(aes(col=`DxOrdinal`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age Acceleration", title='Ordinal Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom') 


pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/dx_vs_age_acceleration_by_region.pdf',height=10,width=10)
full_dx_AA
subset_dx_AA
ordinal_dx_AA
dev.off()

## Age
full_dx_age = ggplot(data=pd, aes(x=Dx,y=Age )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age", title='Full Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 
		
subset_dx_age = ggplot(data=pd[pd$DxOrdinal!="Alz Drop", ], aes(x=Dx,y=Age )) + 
geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age", title='Subset Model') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='bottom') 

ordinal_dx_age = ggplot(data=pd, aes(x=DxOrdinal,y=Age )) + 
geom_point(aes(col=`DxOrdinal`),position = position_jitterdodge(jitter.width=0.4,dodge.width=.85)) +
geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
labs(x='Disease', y = "Age", title='Ordinal Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom') 


pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/dx_vs_age_by_region.pdf',height=10,width=10)
full_dx_age
subset_dx_age
ordinal_dx_age
dev.off()

## aa versus quality
qual_AA1 = ggplot(data=pd, aes(x=negControl_PC1,y=AA_Linear )) + 
geom_point(aes(col=`DxOrdinal`),) +
geom_smooth(se=FALSE,method='lm') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Negative Control PC1', y = "Age Acceleration", title='Full Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom') 
		
qual_AA2 = ggplot(data=pd, aes(x=negControl_PC1,y=AA_Linear )) + 
geom_point(aes(col=`DxOrdinal`),) +
geom_smooth(se=FALSE,method='lm', aes(col=DxOrdinal) ) + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Negative Control PC1', y = "Age Acceleration", title='Full Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom') 
		
qual_AA3 = ggplot(data=pd, aes(x=negControl_PC2,y=AA_Linear )) + 
geom_point(aes(col=`DxOrdinal`),) +
geom_smooth(se=FALSE,method='lm') + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Negative Control PC2', y = "Age Acceleration", title='Full Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom')
		
qual_AA4 = ggplot(data=pd, aes(x=negControl_PC2,y=AA_Linear )) + 
geom_point(aes(col=`DxOrdinal`),) +
geom_smooth(se=FALSE,method='lm', aes(col=DxOrdinal)) + 
facet_wrap(~`Region`,scales='fixed',nrow=1) + 
labs(x='Negative Control PC2', y = "Age Acceleration", title='Full Model') + 
		scale_fill_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		scale_colour_manual(values=c("black",'steelblue',"#ab1323" ) ) + 
		theme(legend.position='bottom') 		
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/qual_age_acceleration_by_region.pdf',height=10,width=10)
qual_AA1
qual_AA2
qual_AA3
qual_AA4
dev.off()
###
library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/exploratory_age_acceleration_heatmap.pdf',height=10,width=10,onefile=T)
pheatmap(cor(pd[,c("AA_Linear",'Age','DNAmAge','negControl_PC1','negControl_PC2','negControl_PC3','negControl_PC4','mMed','xMed',"mean_detection_p",'NeuN_neg','NeuN_pos')]),
		cluster_rows=T, 
		cluster_cols=T,
#		color=col.pal,
		fontsize=14)
		
		pheatmap(abs(cor(pd[,c("AA_Linear",'Age','DNAmAge','negControl_PC1','negControl_PC2','negControl_PC3','negControl_PC4','mMed','xMed',"mean_detection_p",'NeuN_neg','NeuN_pos')])),
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=14)
dev.off()

cor.test(pd$AA_Linear,pd$DNAmAge)
summary(lm(AA_Linear~mMed+uMed+Dx+negControl_PC2+negControl_PC1,data=pd) )

summary(lm(DNAmAge~Age,data=pd) )
summary(lm(DNAmAge~Age+mMed+uMed+negControl_PC2+negControl_PC1+Sex+Dx+Race,data=pd) )
summary(lm(DNAmAge~mMed+uMed+negControl_PC2+negControl_PC1,data=pd) )

cor.test(pd$AA_Linear,pd$mMed)