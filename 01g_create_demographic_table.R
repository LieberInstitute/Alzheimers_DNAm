library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('rdas/cleanSamples_n377_processed_data_postfiltered.rda',verbose=T)

fisher.test(table(pd$APOE4_Dosage[!duplicated(pd$BrNum)&pd$keepList], pd$Dx[!duplicated(pd$BrNum)&pd$keepList]))

fisher.test(table(pd$APOE4_Dosage[pd$Region=="CRB" & pd$keepList], pd$Dx[pd$Region=="CRB" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="DLPFC" & pd$keepList], pd$Dx[pd$Region=="DPLFC" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="HIPPO" & pd$keepList], pd$Dx[pd$Region=="HIPPO" & pd$keepList]))
fisher.test(table(pd$APOE4_Dosage[pd$Region=="ERC" & pd$keepList], pd$Dx[pd$Region=="ERC" & pd$keepList]))

##
library(LIBDpheno)
library(tableone)
pd = cbind(pd, toxicant[[1]][match(brnum(pd$BrNum),toxicant[[1]]$brnumerical),])
pd$APOE4_Dosage = factor(pd$APOE4_Dosage)

pd$DxOrdinal = plyr::revalue(pd$DxOrdinal,c("Alz Drop"="Asymptomatic Alzheimers", "Alz Keep"="Symptomatic Alzheimers") )

t1 <- CreateTableOne(vars = c("Age","Race","Sex",'manner_of_death','cerad','braak','APOE4_Dosage'), strata = c("DxOrdinal","Region"), data = pd)
table1 <- print(t1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
write.csv(table1, file = "csvs/SupplementalTable_summary_of_sample_characteristics_table1.csv")

######
library(ggplot2)
theme_set(theme_bw(base_size=30) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
######
library(tidyverse)
pd$DxOrdinal = plyr::revalue(pd$DxOrdinal,c("Asymptomatic Alzheimers"="Asymptomatic AD", "Symptomatic Alzheimers"="Symptomatic AD") )
pd$DxOrdinal=factor(pd$DxOrdinal,levels=c('Control',"Asymptomatic AD","Symptomatic AD"))
pd=pd[pd$DxOrdinal!="Asymptomatic AD",]

dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(braak)
dat$DxOrdinal = droplevels(dat$DxOrdinal)

a <- ggplot(data=dat, aes(x=braak, y=Age, fill=DxOrdinal)) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_point(aes(col=`DxOrdinal`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 	
		scale_fill_manual(values=c("#ab1323" ) ) +
		scale_colour_manual(values=c("#ab1323" ) ) + 		
		theme(legend.position='none') + labs(x="Braak Stage")
		

dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(cerad)
dat$DxOrdinal = droplevels(dat$DxOrdinal)
		
b <- ggplot(data=dat, aes(x=cerad, y=Age, fill=DxOrdinal)) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_point(aes(col=`DxOrdinal`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 	
		scale_fill_manual(values=c("#ab1323" ) ) +
		scale_colour_manual(values=c("#ab1323" ) ) + 		
		theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = .95, vjust=1))  + 
		labs(x="CERAD Score")

dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(Age)
dat$DxOrdinal = droplevels(dat$DxOrdinal)		
		
c <- ggplot(data=dat, aes(x=DxOrdinal, y=Age,fill=DxOrdinal,col=DxOrdinal) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_jitter(width=0.2)+ theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=1),axis.title.x=element_blank())  + 	
		scale_fill_manual(values=c("black","#ab1323" ) ) +
		scale_colour_manual(values=c("black","#ab1323" ) ) 


dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(brain_weight_gram)
dat$DxOrdinal = droplevels(dat$DxOrdinal)		
		
d <- ggplot(data=dat, aes(x=DxOrdinal, y=brain_weight_gram,fill=DxOrdinal,col=DxOrdinal) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_jitter(width=0.2) + theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=1),axis.title.x=element_blank()) +
		labs(y='Brain Weight (grams)')+ 	
		scale_fill_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ) ) +
		scale_colour_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ) )

dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(braak)
dat$DxOrdinal = droplevels(dat$DxOrdinal)
		
e <- ggplot(data=dat, aes(x=DxOrdinal, y=as.numeric(braak),fill=DxOrdinal,col=DxOrdinal) ) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_jitter(width=0.2)	+ theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=1),axis.title.x=element_blank()) +
		labs(y='Braak Stage')+ 	
		scale_fill_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ) ) +
		scale_colour_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ))

dat = pd[!duplicated(pd$BrNum),]
dat = dat %>% drop_na(cerad)
dat$DxOrdinal = droplevels(dat$DxOrdinal)
		
f <- ggplot(data=dat, aes(x=DxOrdinal, y=as.numeric(cerad),fill=DxOrdinal,col=DxOrdinal )) +
        geom_boxplot(outlier.colour = NA, alpha = 0.5, col='black')  + 
		geom_jitter(width=0.2)	+ theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust=1),axis.title.x=element_blank() ) + 
		labs(y='CERAD Score')+ 	
		scale_fill_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ) ) +
		scale_colour_manual(values=c(Control="black",`Symptomatic AD`="#ab1323" ) )
		
library(gridExtra)
library(grid)
lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))
demo_boxplots <- arrangeGrob(a, b, c,d ,layout_matrix =lay) 
ggsave(demo_boxplots, file="plots/SupplementalFigure_demographic_boxplots.pdf", height=14,width=14,useDingbats=FALSE) 	 	

### Get list of brnums to get CERAD +Braak data for
#pd$braak_or_cerad_missing = ifelse(is.na(pd$braak) | is.na(pd$cerad), "Yes", "No" )
write.csv( pd[!duplicated(pd$BrNum),c('BrNum','braak_or_cerad_missing','source','braak','cerad','DxOrdinal')], file='csvs/Alz_98donors_neuropath_data_available_to_SS.csv',row.names=F )

## t-tests of interest
t.test(Age~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) 
t.test(brain_weight_gram~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) 
t.test(bmi_calculated~DxOrdinal,data=pd[!duplicated(pd$BrNum) & pd$DxOrdinal!="Alz Drop",]) #
			