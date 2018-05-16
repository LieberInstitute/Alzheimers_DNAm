### Check overlap of genetic risk with 
ad_genetic_risk_loci=read.csv('/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_loci_cpgs.csv') 

### 
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')

cpg_in_loci = cbind(allStats[as.character(ad_genetic_risk_loci$Name),],ad_genetic_risk_loci[,-1] )
table(is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val))
cpg_in_loci = cpg_in_loci[!is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val),]
## Create summary table
library(dplyr)
library(tidyr)
cpg_in_loci = group_by(cpg_in_loci,nearestGeneToSNP )
genetic_risk_DMP_summary = summarize(cpg_in_loci, 
									 GeneticMarker = unique(snpRsNum),
									 Chromosome = unique(chr_hg38), 
									 Covered_Region = paste0(min(pos_hg38),":",max(pos_hg38)) , 
									 Region_Length = max(pos_hg38) - min(pos_hg38), 
									 Number_Of_CpGs = n(), 
									 n_mainEffect_nomP05= sum(Primary_subset_mainEffect_P.Value<0.05),
									 n_mainEffect_FDR05= sum(Primary_subset_mainEffect_adj.P.Val<0.05),
									 n_interaction_nomP05=sum(Primary_subset_interactionEffect_P.Value<0.05),
									 n_interactionEffect_FDR05= sum(Primary_subset_interactionEffect_adj.P.Val<0.05)) %>%
									 as.data.frame()
write.csv(genetic_risk_DMP_summary, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/genetic_risk_dmp_overlap_summary.csv',row.names=F )									 

### Make plots
library(ggplot2)
cpg_in_loci_df = as.data.frame(cpg_in_loci)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
				 
pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/AD_gwas_risk_loci_methylation_signal_mainEffect.pdf')
for (i in levels(cpg_in_loci_df$nearestGeneToSNP) ){
dat= cpg_in_loci_df[cpg_in_loci_df$nearestGeneToSNP%in%i, ]

dat$sig = factor(ifelse(dat$Primary_subset_mainEffect_adj.P.Val<0.05, "FDR<0.05","FDR>0.05"))
levels(dat$sig) <- c("FDR<0.05","FDR>0.05")
colScale <- scale_colour_manual(name = "sig",values = c(`FDR>0.05`='black',`FDR>0.05`='red') )

a = ggplot(data=dat,aes(x=pos_hg38/1e3, y=-log10(Primary_subset_mainEffect_P.Value) ) ) + geom_point() + labs(x=paste0('CpG position on ',unique(dat$chr_hg38),' in kilobases (hg38)'),y='-log10(mainEffect MethP)',title=paste0('Differential methylation at ',as.character(i) ,' locus') ) + scale_y_continuous(limits = c(0, max(-log10(cpg_in_loci_df$Primary_subset_mainEffect_P))+0.5 ), expand = c(0, 0)) 
print(a)
}
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/Alz/plots/AD_gwas_risk_loci_methylation_signal_interactionEffect.pdf')
for (i in levels(cpg_in_loci_df$nearestGeneToSNP) ){
dat= cpg_in_loci_df[cpg_in_loci_df$nearestGeneToSNP%in%i, ]

dat$sig = factor(ifelse(dat$Primary_subset_interactionEffect_adj.P.Val<0.05, "FDR<0.05","FDR>0.05"))
levels(dat$sig) <- c("FDR>0.05","FDR<0.05")

a = ggplot(data=dat,aes(x=pos_hg38/1e3, y=-log10(Primary_subset_interactionEffect_P.Value) ) ) + geom_point() + labs(x=paste0('CpG position on ',unique(dat$chr_hg38),' in kilobases (hg38)'),y='-log10(interactionEffect MethP)',title=paste0('Differential methylation at ',as.character(i) ,' locus') ) + scale_y_continuous(limits = c(0, max(-log10(cpg_in_loci_df$Primary_subset_interactionEffect_P))+0.5 ), expand = c(0, 0)) 
print(a)
}
dev.off()

### Check enrichment
allStats$within_genetic_risk_loci = allStats$Name%in% ad_genetic_risk_loci$Name

fisher.test(table(riskLoci=allStats$within_genetic_risk_loci,DMP=allStats$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_mainEffect_P.Value<0.05 ) )

fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_P.Value<0.05 ) )

################################# Re-run all analysis with LD block loci
ad_genetic_risk_loci = read.csv('/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_LD_block_loci_cpgs.csv')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')
cpg_in_loci = cbind(allStats[as.character(ad_genetic_risk_loci$Name),],ad_genetic_risk_loci[,-1] )
table(is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val))
cpg_in_loci = cpg_in_loci[!is.na(cpg_in_loci$Primary_subset_mainEffect_adj.P.Val),]

## Create summary table
library(dplyr)
library(tidyr)
cpg_in_loci = group_by(cpg_in_loci,indexSnp_nearestGene )
genetic_risk_DMP_summary = summarize(cpg_in_loci, 
									 GeneticMarker = unique(indexSnp_rsNum),
									 Chromosome = unique(chr_hg38), 
									 Covered_Region = paste0(min(pos_hg38),":",max(pos_hg38)) , 
									 Region_Length = max(pos_hg38) - min(pos_hg38), 
									 Number_Of_CpGs = n(), 
									 n_mainEffect_nomP05= sum(Primary_subset_mainEffect_P.Value<0.05),
									 n_mainEffect_FDR05= sum(Primary_subset_mainEffect_adj.P.Val<0.05),
									 n_interaction_nomP05=sum(Primary_subset_interactionEffect_P.Value<0.05),
									 n_interactionEffect_FDR05= sum(Primary_subset_interactionEffect_adj.P.Val<0.05)) %>%
									 as.data.frame()
write.csv(genetic_risk_DMP_summary, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/genetic_risk_dmp_overlap_summary_ld_block_loci.csv',row.names=F )

### Check enrichment
allStats$within_genetic_risk_loci = allStats$Name%in% ad_genetic_risk_loci$Name
######## Main effect
### All
fisher.test(table(riskLoci=allStats$within_genetic_risk_loci,DMP=allStats$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_mainEffect_P.Value<0.05 ) )
### hypermethylated
hyper = allStats[allStats$Primary_subset_mainEffect_logFC>0,]
fisher.test(table(riskLoci=hyper$within_genetic_risk_loci,DMP=hyper$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(hyper$within_genetic_risk_loci,hyper$Primary_subset_mainEffect_P.Value<0.05 ) )

### hypomethylated
hypo = allStats[allStats$Primary_subset_mainEffect_logFC<0,]
fisher.test(table(riskLoci=hypo$within_genetic_risk_loci,DMP=hypo$Primary_subset_mainEffect_adj.P.Val<0.05 ) )
fisher.test(table(hypo$within_genetic_risk_loci,hypo$Primary_subset_mainEffect_P.Value<0.05 ) )

######## interaction effect
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_adj.P.Val<0.05 ) )
fisher.test(table(allStats$within_genetic_risk_loci,allStats$Primary_subset_interactionEffect_P.Value<0.05 ) )

