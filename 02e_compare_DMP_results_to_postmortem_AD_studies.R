library(RColorBrewer)
library(pheatmap)
col.pal = brewer.pal(9,"Blues")
setwd('/dcl01/lieber/ajaffe/Steve/Alz')

#load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')
###

#Lunnon_CROSS_DMPs= openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/nn.3782-S7.xlsx',startRow=3)

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

###
###
jager_cpgs = intersect(Jager_S2[,"Target ID"], mergedStats$Name)
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
##### Replication stats by all region		
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_allRegion.rda')
										 
allRegion_DMP_Replication=group_by(allRegion_tidyStats,Interaction, Sensitivity, Model) %>% summarise(  Jager_Replicate_N= check_replication(
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

write.csv(allRegion_DMP_Replication, '/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_DMP_Results_Replication.csv',row.names=F)

##### Replication stats by all region		
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/tidyStats_caseControl_DMC_singleRegion.rda')
										 
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
write.csv(singleRegion_DMP_Replication, '/dcl01/lieber/ajaffe/Steve/Alz/csvs/singleRegion_DMP_Results_Replication.csv',row.names=F)

#### Check Watson
Watson=openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/Alz/13073_2015_258_MOESM6_ESM.xlsx',start=3)
## Check the jager cpgs in our data
shared_cpgs = intersect(Watson[,"ProbeID"], mergedStats$Name)
watson_cpgs = mergedStats[mergedStats$Name%in%shared_cpgs,]


watson_matrix = -log10(watson_cpgs[,grepl("P.Value",colnames(watson_cpgs) )] )
pdf("plots/pheatmap_watson_cpg_pvalues_across_models.pdf",h=30,w=30,onefile=TRUE)
pheatmap(watson_matrix, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)
		
pheatmap(watson_matrix[,-grep("ALL",colnames(watson_matrix) )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)	
		
pheatmap(watson_matrix[,grep("ALL",colnames(watson_matrix) )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)		
dev.off()

## Check Yu

yu=tabulizer::extract_tables('/dcl01/lieber/ajaffe/Steve/Alz/NOI140084supp1_prod.pdf')
##
## Check the Lunnon cpgs in our data
Lunnon_data= getGEO("GSE59685")
Lunnon_pd = pData(Lunnon_data[[1]])
filePaths = getGEOSuppFiles("GSE59685")
untar(rownames(filePaths)[1], exdir='GSE59685')
Lunnon_beta=read.table(gzfile("/dcl01/lieber/ajaffe/Steve/Alz/GSE59685/GSE59685_betas.csv.gz"),sep=",")  


## Check the jager cpgs in our data
shared_cpgs = intersect(Jager_S2[,"Target ID"], mergedStats$Name)
jager_cpgs = mergedStats[mergedStats$Name%in%shared_cpgs,]

## Analyze the jager p values
jager_matrix_p = jager_cpgs[,grepl("P.Value",colnames(jager_cpgs) )]
jager_matrix_p = cbind(jager_matrix_p, Jager_S2[match(rownames(jager_cpgs), Jager_S2[,"Target ID"]),grep("p-value",colnames(Jager_S2))])
sapply(jager_matrix_p, function(x) (sum(x<0.05)/70) *100 )
jager_matrix_pLog = -log10(jager_matrix_p)

pdf("plots/pheatmap_jager_cpg_pvalues_across_models.pdf",h=30,w=30,onefile=TRUE)
pheatmap(jager_matrix_pLog, 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)
		
pheatmap(jager_matrix_pLog[,-grep("ALL",colnames(jager_matrix_pLog) )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)	
		
pheatmap(jager_matrix_pLog[,grep("ALL",colnames(jager_matrix_pLog) )], 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)		
dev.off()

pdf("plots/pheatmap_jager_cpg_pvalue_model_correlation.pdf",h=30,w=30,onefile=TRUE)
pheatmap(cor(jager_matrix_pLog), 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)	
dev.off()

## Analyze the jager model estimates
jager_matrix_est = jager_cpgs[,grepl("logFC",colnames(jager_cpgs) )]
jager_matrix_est = cbind(jager_matrix_est, Jager_S2[match(rownames(jager_cpgs), Jager_S2[,"Target ID"]),grep("Est",colnames(Jager_S2))])

jager_matrix_est_rescaled=jager_matrix_est
jager_matrix_est_rescaled[11:14] = jager_matrix_est[11:14]/100
pdf("plots/pheatmap_jager_cpg_estimates_across_models.pdf",h=30,w=30,onefile=TRUE)
pheatmap(jager_matrix_est_rescaled, 
		cluster_rows=T, 
		cluster_cols=T,
	#	color=col.pal,
		fontsize=20)
pheatmap(cor(jager_matrix_est_rescaled), 
		cluster_rows=T, 
		cluster_cols=T,
		color=col.pal,
		fontsize=20)
dev.off()


## Check the jager cpgs in our data
