###### Test all-region models
library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')
##
testOntologies = function(sigCpg, allCpG) {
cat(".")
out <- tryCatch ( {
gstGO <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=allCpG, collection="GO")
gstKEGG <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=allCpG, collection="KEGG")
return( list(GO=gstGO,KEGG=gstKEGG) ) },
error=function(cond) {
            # Choose a return value in case of error
            return(NULL)
        }
)
return(out) }

############# FDR 5% split ##########
tests = grep("adj.P.Val", colnames(allStats),value=TRUE )
mainTests = grep("Primary_subset_main",tests, value=T)

colSums(allStats[,mainTests,drop=F] <0.05)
##
sigStatList = lapply(mainTests, function(test) { 
  x = allStats[which(allStats[,test]<0.05), c("Name",gsub("_adj.P.Val","_logFC",test))]
  x$Test = gsub("_adj.P.Val","", test)
  colnames(x) <-c("Name","LFC","Test")
  return(x)} )
  
sigStats = do.call("rbind", sigStatList)
sigStats$Sign = sign(sigStats$LFC)
sigStats$lab =paste0(sigStats$Sign, "_", sigStats$Test)

sigCpg = split(sigStats$Name, sigStats$lab)
allCpG=unique(allStats$Name)

##

allRegion_GeneOntology = lapply(sigCpg, testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","", names(sigCpg) )
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
## Reordering
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder( jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1)) ]
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(gsub(".*1_","", names(allRegion_GeneOntology_list_split) ) )]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_main_GO_Analysis_FDR05_DMC_case_control_stats_Split.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_main_GO_Analysis_FDR05_DMP_Results_Split.xlsx')

openxlsx::write.xlsx(allRegion_GeneOntology_list_split[grep("-1_", names(allRegion_GeneOntology_list_split))], file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_main_GO_Analysis_FDR05_DMP_Results_hypomethylated.xlsx')
openxlsx::write.xlsx(allRegion_GeneOntology_list_split[-grep("-1_", names(allRegion_GeneOntology_list_split))], file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_main_GO_Analysis_FDR05_DMP_Results_hypermethylated.xlsx')

############# FDR 5% interaction ##########
interactionTests = "Primary_subset_interactionEffect_adj.P.Val"
colSums(allStats[,interactionTests,drop=F] <0.05)
allCpG=unique(allStats$Name)
sigCpg= lapply(interactionTests, function(test) {allStats[allStats[,test]<.05, 'Name']})

interaction_GeneOntology = lapply(sigCpg,testOntologies, allCpG)
names(interaction_GeneOntology) <- "Interaction"

##
interaction_GeneOntology_list= unlist(interaction_GeneOntology,recursive=F)
interaction_GeneOntology_list_split = c( interaction_GeneOntology_list[grep("KEGG",names(interaction_GeneOntology_list) )],unlist(lapply(interaction_GeneOntology_list[grep("GO",names(interaction_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )

interaction_GeneOntology_list_split= lapply(interaction_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
interaction_GeneOntology_list_split= lapply(interaction_GeneOntology_list_split, tibble::rownames_to_column  )
##save
openxlsx::write.xlsx(interaction_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/interaction_GO_Analysis_FDR05_DMP_Results_dupCor.xlsx')



########################################### Everything else


##### FDR 10 full ####################
allCpG=unique(allStats$Name)
sigCpg= lapply(tests, function(test) {allStats[allStats[,test]<.10, 'Name']})

allRegion_GeneOntology = lapply(sigCpg,testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","",tests)
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("interaction","int",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1))]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_GO_Analysis_FDR10_DMC_case_control_stats_Unsplit.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_GO_Analysis_FDR10_DMP_Results_Unsplit.xlsx')

##### FDR 5 full ####################
allCpG=unique(allStats$Name)
sigCpg= lapply(tests, function(test) {allStats[allStats[,test]<.05, 'Name']})

allRegion_GeneOntology = lapply(sigCpg,testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","",tests)
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("interaction","int",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1))]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_GO_Analysis_FDR10_DMC_case_control_stats_Unsplit.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_GO_Analysis_FDR10_DMP_Results_Unsplit.xlsx')


##### FDR 1 full ####################
allCpG=unique(allStats$Name)
sigCpg= lapply(tests, function(test) {allStats[allStats[,test]<.01, 'Name']})

allRegion_GeneOntology = lapply(sigCpg,testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","",tests)
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("interaction","int",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1))]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_GO_Analysis_FDR1_DMC_case_control_stats_Unsplit.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_GO_Analysis_FDR1_DMP_Results_Unsplit.xlsx')

############# FDR 10 split ##########
tests = grep("adj.P.Val", colnames(allStats),value=TRUE )
mainTests = grep("main",tests, value=T)

colSums(allStats[,mainTests] <0.10)
##
sigStatList = lapply(mainTests, function(test) { 
  x = allStats[which(allStats[,test]<0.1), c("Name",gsub("_adj.P.Val","_logFC",test))]
  x$Test = gsub("_adj.P.Val","", test)
  colnames(x) <-c("Name","LFC","Test")
  return(x)} )
  
sigStats = do.call("rbind", sigStatList)
sigStats$Sign = sign(sigStats$LFC)
sigStats$lab =paste0(sigStats$Sign, "_", sigStats$Test)

sigCpg = split(sigStats$Name, sigStats$lab)
allCpG=unique(allStats$Name)

##

allRegion_GeneOntology = lapply(sigCpg, testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","", names(sigCpg) )
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
## Reordering
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder( jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1)) ]
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(gsub(".*1_","", names(allRegion_GeneOntology_list_split) ) )]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_GO_Analysis_FDR10_DMC_case_control_stats_Split.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_GO_Analysis_FDR10_DMP_Results_Split.xlsx')

############# FDR 1 split ##########
tests = grep("adj.P.Val", colnames(allStats),value=TRUE )
mainTests = grep("main",tests, value=T)

colSums(allStats[,mainTests] <0.01)
##
sigStatList = lapply(mainTests, function(test) { 
  x = allStats[which(allStats[,test]<0.01), c("Name",gsub("_adj.P.Val","_logFC",test))]
  x$Test = gsub("_adj.P.Val","", test)
  colnames(x) <-c("Name","LFC","Test")
  return(x)} )
  
sigStats = do.call("rbind", sigStatList)
sigStats$Sign = sign(sigStats$LFC)
sigStats$lab =paste0(sigStats$Sign, "_", sigStats$Test)

sigCpg = split(sigStats$Name, sigStats$lab)
allCpG=unique(allStats$Name)

##

allRegion_GeneOntology = lapply(sigCpg, testOntologies, allCpG)
names(allRegion_GeneOntology) <- gsub("_adj.P.Val","", names(sigCpg) )
names(allRegion_GeneOntology) = gsub("Effect","",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Primary","Unadj",names(allRegion_GeneOntology))
names(allRegion_GeneOntology) = gsub("Sensitivity","Adj",names(allRegion_GeneOntology))

##
allRegion_GeneOntology_list= unlist(allRegion_GeneOntology,recursive=F)
allRegion_GeneOntology_list_split = c( allRegion_GeneOntology_list[grep("KEGG",names(allRegion_GeneOntology_list) )],unlist(lapply(allRegion_GeneOntology_list[grep("GO",names(allRegion_GeneOntology_list) )], function(x) split(x, x$Ont)),recursive=F) )
## Reordering
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder( jaffelab::ss(names(allRegion_GeneOntology_list_split), "\\.",1)) ]
allRegion_GeneOntology_list_split = allRegion_GeneOntology_list_split[gtools::mixedorder(gsub(".*1_","", names(allRegion_GeneOntology_list_split) ) )]


allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, function(x) x[order(x$'FDR',x$'P.DE'),] )
allRegion_GeneOntology_list_split= lapply(allRegion_GeneOntology_list_split, tibble::rownames_to_column  )
##save
save(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_GO_Analysis_FDR1_DMC_case_control_stats_Split.rda')	
openxlsx::write.xlsx(allRegion_GeneOntology_list_split, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/allRegion_GO_Analysis_FDR1_DMP_Results_Split.xlsx')

### nominal pathway analysis
tests = grep("adj.P.Value", colnames(mergedStats),value=TRUE )
tests = tests[!grepl("full",tests)]

##1e-3
colSums(mergedStats[,tests] <0.001)



GO_Analysis = lapply(tests, testOntologies)
names(GO_Analysis) <- gsub("_P.Value","",tests)
GO_Analysis_ss = GO_Analysis[lengths(GO_Analysis)==2]
GO_Analysis_ss= unlist(GO_Analysis_ss,recursive=F)
GO_Analysis_ss= lapply(GO_Analysis_ss, function(x) x[order(x$'FDR',x$'P.DE'),] )

save(GO_Analysis_ss, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/GO_Analysis_p1e3_DMC_case_control_stats.rda')	
openxlsx::write.xlsx(GO_Analysis_ss, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/GO_Analysis_p1e3_DMP_Results.xlsx')

##1e-4
colSums(mergedStats[,tests] <1e-4)

GO_Analysis = lapply(tests, testOntologies,thresh=1e-4)
names(GO_Analysis) <- gsub("_P.Value","",tests)
GO_Analysis_ss = GO_Analysis[lengths(GO_Analysis)==2]
GO_Analysis_ss= unlist(GO_Analysis_ss,recursive=F)
GO_Analysis_ss= lapply(GO_Analysis_ss, function(x) x[order(x$'FDR',x$'P.DE'),] )

save(GO_Analysis_ss, file='/dcl01/lieber/ajaffe/Steve/Alz/rdas/GO_Analysis_p1e4_DMC_case_control_stats.rda')	
openxlsx::write.xlsx(GO_Analysis_ss, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/GO_Analysis_p1e4_DMP_Results.xlsx')


######### Using adjusted pvalue cutoff
tests = grep("adj.P.Val", colnames(mergedStats),value=TRUE )
testOntologies = function(tests) {
cat(".")
sigCpg= mergedStats[mergedStats[,tests]<0.10, 'Name']
out <- tryCatch ( {
gstGO <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=mergedStats[,'Name'], collection="GO")
gstKEGG <- missMethyl::gometh(sig.cpg=sigCpg, all.cpg=mergedStats[,'Name'], collection="KEGG")
return( list(GO=gstGO,KEGG=gstKEGG) ) },
error=function(cond) {
            # Choose a return value in case of error
            return(NULL)
        }
)
return(out) }

GO_Analysis = lapply(tests, testOntologies)
names(GO_Analysis) <- gsub("_adj.P.Val","",tests)
GO_Analysis_ss = GO_Analysis[lengths(GO_Analysis)==2]
GO_Analysis_ss= unlist(GO_Analysis_ss,recursive=F)
save(GO_Analysis_ss, '/dcl01/lieber/ajaffe/Steve/Alz/rdas/GO_Analysis_DMC_case_control_stats.rda')	
openxlsx::write.xlsx(GO_Analysis_ss, file='/dcl01/lieber/ajaffe/Steve/Alz/csvs/GO_Analysis_DMP_Results.xlsx')
