setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
load('/dcl01/ajaffe/data/lab/libd_alzheimers/grant_analysis_hg38/LIBD_AD_pd_hg38.Rdata')

pd$Region = plyr::revalue(pd$Region, c("CB"="CRB") )

### Determine case-control statistics
p
summary( aov((pd[,'mitoRate']~pd[,'Region'] ) ))
summary( aov((pd[,'totalAssignedGene']~pd[,'Region'] ) ))
summary( aov((pd[,'RIN']~pd[,'Region'] ) ))

#Test for mitochondrial map rate associations
t.test(pd[pd$Region=="CRB",'mitoRate']~pd[pd$Region=="CRB",'Dx'] )
t.test(pd[pd$Region=="ERC",'mitoRate']~pd[pd$Region=="ERC",'Dx'] )
t.test(pd[pd$Region=="DLPFC",'mitoRate']~pd[pd$Region=="DLPFC",'Dx'] )
t.test(pd[pd$Region=="HIPPO",'mitoRate']~pd[pd$Region=="HIPPO",'Dx'] )

#Test for GAR associations
t.test(pd[pd$Region=="CRB",'totalAssignedGene']~pd[pd$Region=="CRB",'Dx'] )
t.test(pd[pd$Region=="ERC",'totalAssignedGene']~pd[pd$Region=="ERC",'Dx'] )
t.test(pd[pd$Region=="DLPFC",'totalAssignedGene']~pd[pd$Region=="DLPFC",'Dx'] )
t.test(pd[pd$Region=="HIPPO",'totalAssignedGene']~pd[pd$Region=="HIPPO",'Dx'] )

#Test for RIN associations
t.test(pd[pd$Region=="CRB",'RIN']~pd[pd$Region=="CRB",'Dx'] )
t.test(pd[pd$Region=="ERC",'RIN']~pd[pd$Region=="ERC",'Dx'] )
t.test(pd[pd$Region=="DLPFC",'RIN']~pd[pd$Region=="DLPFC",'Dx'] )
t.test(pd[pd$Region=="HIPPO",'RIN']~pd[pd$Region=="HIPPO",'Dx'] )

### Produce case-control boxplots
library(RColorBrewer)
library(ggplot2)
theme_set(theme_bw(base_size=40) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

a = ggplot(pd, aes_string(x = 'Region', y = 'totalAssignedGene', fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y='GAR',x='Diagnosis') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30,colour='black'))
		
b = ggplot(pd, aes_string(x = 'Region', y = 'mitoRate', fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y='mitoMap Rate',x='Diagnosis') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30,colour='black'))

c = ggplot(pd, aes_string(x = 'Region', y = 'RIN', fill='Dx')) +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_point(aes(col=`Dx`),position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
		labs(y='RIN',x='Diagnosis') + 
		scale_fill_manual(values=c("black","#ab1323" ) ) + 
		scale_colour_manual(values=c("black","#ab1323" ) ) + 
		theme(legend.position='none') +
		theme(axis.title.x=element_blank(), axis.text.x=element_text(size=30,colour='black'))		

library(gridExtra)
library(grid)
options(bitmapType="cairo")
X11(type = "cairo")
qcRNAseq = grid.arrange(a, b, c, nrow=3)

ggsave(qcRNAseq, file="qc/RNA_seq_metrics_caseControl_byRegion_boxplots.pdf", height=16,width=12,useDingbats=FALSE) 	 	

		
pdf('qc/RNA_seq_metrics_caseControl_byRegion_boxplots.pdf',height=16,width=12)
a
b
c
dev.off()