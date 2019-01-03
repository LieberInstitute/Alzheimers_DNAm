library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')

###  Check Jager genes
library(docxtractr)
Jager_et_al <- read_docx("/dcl01/lieber/ajaffe/Steve/Alz/NIHMS614693-supplement-supp_tables_only.docx")
Jager_et_al_tables = docx_extract_all_tbls(Jager_et_al)
Jager_S2 = Jager_et_al_tables[[3]]
Jager_S2_colnames =  c( colnames(Jager_S2)[1:3], "NP burden Model 1 Est","NP burden Model 1 p-value", "NP burden Model 3 Est","NP burden Model 3 p-value",	"AD Est", "AD p-value",	"Braak score (PFC) Est", "Braak score (PFC) p-value", "Genes within 50 kb")
Jager_S2 = Jager_S2[-(1:2), ]
colnames(Jager_S2) = Jager_S2_colnames
Jager_S2[,4:11] = sapply(Jager_S2[,4:11], as.numeric)

#### ERC #####
Lunnon_ERC_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S4.xlsx',startRow=2)
prefixes=rep(colnames(Lunnon_ERC_DMPs)[-grep("X",colnames(Lunnon_ERC_DMPs))], each=2)
column_names = Lunnon_ERC_DMPs[1,-1]
Lunnon_ERC_DMPs =Lunnon_ERC_DMPs[-1,-1]
colnames(Lunnon_ERC_DMPs) <- column_names
colnames(Lunnon_ERC_DMPs) = gsub("P value ","P", colnames(Lunnon_ERC_DMPs))
colnames(Lunnon_ERC_DMPs) = gsub(colnames(Lunnon_ERC_DMPs)[7],"estimate", colnames(Lunnon_ERC_DMPs))
k=min( grep("P$", colnames(Lunnon_ERC_DMPs)) ):ncol(Lunnon_ERC_DMPs)
colnames(Lunnon_ERC_DMPs)[k] <- paste0( prefixes,"_",colnames(Lunnon_ERC_DMPs)[k])
Lunnon_ERC_DMPs[,k]=sapply(Lunnon_ERC_DMPs[,k],as.numeric)

#### STG #####
Lunnon_STG_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S5.xlsx',startRow=2)
prefixes=rep(colnames(Lunnon_STG_DMPs)[-grep("X",colnames(Lunnon_STG_DMPs))], each=2)
column_names = Lunnon_STG_DMPs[1,-1]
Lunnon_STG_DMPs =Lunnon_STG_DMPs[-1,-1]
colnames(Lunnon_STG_DMPs) <- column_names
colnames(Lunnon_STG_DMPs) = gsub("P value ","P", colnames(Lunnon_STG_DMPs))
colnames(Lunnon_STG_DMPs) = gsub(colnames(Lunnon_STG_DMPs)[7],"estimate", colnames(Lunnon_STG_DMPs))
k=min( grep("P$", colnames(Lunnon_STG_DMPs)) ):ncol(Lunnon_STG_DMPs)
colnames(Lunnon_STG_DMPs)[k] <- paste0( prefixes,"_",colnames(Lunnon_STG_DMPs)[k])
Lunnon_STG_DMPs[,k]=sapply(Lunnon_STG_DMPs[,k],as.numeric)

#### PFC #####
Lunnon_PFC_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S6.xlsx',startRow=2)
prefixes=rep(colnames(Lunnon_PFC_DMPs)[-grep("X",colnames(Lunnon_PFC_DMPs))], each=2)
column_names = Lunnon_PFC_DMPs[1,-1]
Lunnon_PFC_DMPs =Lunnon_PFC_DMPs[-1,-1]
colnames(Lunnon_PFC_DMPs) <- column_names
colnames(Lunnon_PFC_DMPs) = gsub("P value ","P", colnames(Lunnon_PFC_DMPs))
colnames(Lunnon_PFC_DMPs) = gsub(colnames(Lunnon_PFC_DMPs)[7],"estimate", colnames(Lunnon_PFC_DMPs))
k=min( grep("P$", colnames(Lunnon_PFC_DMPs)) ):ncol(Lunnon_PFC_DMPs)
colnames(Lunnon_PFC_DMPs)[k] <- paste0( prefixes,"_",colnames(Lunnon_PFC_DMPs)[k])
Lunnon_PFC_DMPs[,k]=sapply(Lunnon_PFC_DMPs[,k],as.numeric)

#### CRB #####
Lunnon_CRB_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S7.xlsx',startRow=2)
prefixes=rep(colnames(Lunnon_CRB_DMPs)[-grep("X",colnames(Lunnon_CRB_DMPs))], each=2)
column_names = Lunnon_CRB_DMPs[1,-1]
Lunnon_CRB_DMPs =Lunnon_CRB_DMPs[-1,-1]
colnames(Lunnon_CRB_DMPs) <- column_names
colnames(Lunnon_CRB_DMPs) = gsub("P value ","P", colnames(Lunnon_CRB_DMPs))
colnames(Lunnon_CRB_DMPs) = gsub(colnames(Lunnon_CRB_DMPs)[7],"estimate", colnames(Lunnon_CRB_DMPs))
k=min( grep("P$", colnames(Lunnon_CRB_DMPs)) ):ncol(Lunnon_CRB_DMPs)
colnames(Lunnon_CRB_DMPs)[k] <- paste0( prefixes,"_",colnames(Lunnon_CRB_DMPs)[k])
Lunnon_CRB_DMPs[,k]=sapply(Lunnon_CRB_DMPs[,k],as.numeric)

#### Cross-region #####
Lunnon_crossCortex_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S8.xlsx',startRow=2)
column_names = Lunnon_crossCortex_DMPs[1,]
Lunnon_crossCortex_DMPs =Lunnon_crossCortex_DMPs[-1,]
colnames(Lunnon_crossCortex_DMPs) <- column_names
colnames(Lunnon_crossCortex_DMPs) = gsub("P value ","P", colnames(Lunnon_crossCortex_DMPs))
Lunnon_crossCortex_DMPs = Lunnon_crossCortex_DMPs[,1:7] #drop cols we don't care about
Lunnon_crossCortex_DMPs[,6:7]=sapply(Lunnon_crossCortex_DMPs[,6:7],as.numeric)
colnames(Lunnon_crossCortex_DMPs)[6:7] <- c('Fisher P', 'Brown P')
##### Replication stats by all region		
load('rdas/tidyStats_caseControl_DMC_allRegion.rda')

###
###
jager_cpgs = intersect(Jager_S2[,"Target ID"], unique(allRegion_tidyStats$Name) )
rownames(Jager_S2) = Jager_S2$`Target ID`

check_replication= function(discovery_names = a1,
							discovery_effect_size= a2,
							replication_names = a3,
							replication_p = a4,
							replication_effect_size = a5) {
shared_names = unique(intersect(discovery_names, replication_names))
disc_i=match(shared_names, discovery_names)
rep_i=match(shared_names, replication_names)
							
N = sum(replication_p[rep_i] <0.05 & sign(discovery_effect_size[disc_i]) == sign(replication_effect_size[rep_i]),na.rm=F )

percent = N*100/length(shared_names)
return( c(N,percent) )
}
myMin = function(x) {min(x)+1}
####
										 
allRegion_DMP_Replication=group_by(allRegion_tidyStats,Interaction, CellTypeAdjustment, Model) %>% summarise(  Jager_Replicate_N= check_replication(
												 discovery_names=Jager_S2[,"Target ID"],
												 discovery_effect_size=Jager_S2[,"AD Est"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Jager_Replicate_Percent = check_replication(
												 discovery_names=Jager_S2[,"Target ID"],
												 discovery_effect_size=Jager_S2[,"AD Est"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2], 
												 
												 Lunnon_ERC_Replicate_N = check_replication(
												 discovery_names=Lunnon_ERC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_ERC_DMPs[,"EC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_ERC_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_ERC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_ERC_DMPs[,"EC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_STG_Replicate_N = check_replication(
												 discovery_names=Lunnon_STG_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_STG_DMPs[,"STG_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_STG_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_STG_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_STG_DMPs[,"STG_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_PFC_Replicate_N = check_replication(
												 discovery_names=Lunnon_PFC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_PFC_DMPs[,"PFC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_PFC_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_PFC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_PFC_DMPs[,"PFC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_CRB_Replicate_N = check_replication(
												 discovery_names=Lunnon_CRB_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_CRB_DMPs[,"CER_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_CRB_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_CRB_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_CRB_DMPs[,"CER_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2]) %>% as.data.frame()

write.csv(allRegion_DMP_Replication, 'csvs/SupplementalTable_allRegion_DMP_Results_Replication.csv',row.names=F)

##### Replication stats by all region		
load('rdas/tidyStats_caseControl_DMC_singleRegion.rda')
										 
singleRegion_DMP_Replication=group_by(singleRegion_tidyStats,Region, Model, CellTypeAdjustment) %>% summarise(  Jager_Replicate_N= check_replication(
												 discovery_names=Jager_S2[,"Target ID"],
												 discovery_effect_size=Jager_S2[,"AD Est"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Jager_Replicate_Percent = check_replication(
												 discovery_names=Jager_S2[,"Target ID"],
												 discovery_effect_size=Jager_S2[,"AD Est"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2], 
												 
												 Lunnon_ERC_Replicate_N = check_replication(
												 discovery_names=Lunnon_ERC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_ERC_DMPs[,"EC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_ERC_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_ERC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_ERC_DMPs[,"EC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_STG_Replicate_N = check_replication(
												 discovery_names=Lunnon_STG_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_STG_DMPs[,"STG_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_STG_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_STG_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_STG_DMPs[,"STG_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_PFC_Replicate_N = check_replication(
												 discovery_names=Lunnon_PFC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_PFC_DMPs[,"PFC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_PFC_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_PFC_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_PFC_DMPs[,"PFC_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2],

												 Lunnon_CRB_Replicate_N = check_replication(
												 discovery_names=Lunnon_CRB_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_CRB_DMPs[,"CER_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[1],
												 Lunnon_CRB_Replicate_Percent = check_replication(
												 discovery_names=Lunnon_CRB_DMPs[,"Probe"],
												 discovery_effect_size=Lunnon_CRB_DMPs[,"CER_estimate"],
												 replication_names=Name,
												 replication_p=P.Value,
												 replication_effect_size=logFC)[2]) %>% as.data.frame()
write.csv(singleRegion_DMP_Replication, 'csvs/SupplementalTable_singleRegion_DMP_Results_Replication.csv',row.names=F)

################ Making venn diagram at CpG level
load('rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')
ourSigCpG = allRegion_mergedStats[allRegion_mergedStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05,'Name']
LunnonCrossCortextCpG = Lunnon_crossCortex_DMPs[,'Probe']
JagerCpG = Jager_S2[,'Target ID']

shared = Reduce(intersect, TopHits)
allRegion_mergedStats[shared,c('Name','within10kb_geneSymbol_gencode_hg38.1','ALL_subset_mainEffect_P.Value')]
a
TopHits = list(LIBD=ourSigCpG, Lunnon_et_al=LunnonCrossCortextCpG, DeJager_et_al=JagerCpG)
####
setwd('/dcl01/lieber/ajaffe/Steve/Alz/Paper')
library(VennDiagram)
## sig
venn_overlap_cpgs = venn.diagram(x = TopHits,
							category.names = c('LIBD','Lunnon et al', 'De Jager et al'),
							filename = NULL,
							fill = c('red', 'blue','green'),
							#cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
							cex=5, cat.cex=2
							)							
grid.draw(venn_overlap_cpgs)

pdf(file='plots/SupplementalFigure_overlap_between_reported_CpGs.pdf',height=12,width=12)
    grid.draw(venn_overlap_cpgs)
dev.off()
