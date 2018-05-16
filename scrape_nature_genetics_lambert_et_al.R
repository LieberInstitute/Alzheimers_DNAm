## Scrape table from nature genetics article
library("rvest")
url <- "https://www.nature.com/articles/ng.2802/tables/2"
risk_loci <- url %>%
  read_html() %>%
  html_nodes(xpath='//*[@id="content"]/div/div/figure/div[1]/div/div[1]/table') %>%
  html_table(fill=T)
risk_loci = risk_loci[[1]]

v=which(apply(risk_loci,1, function(x) length(unique(unlist(x))) )==1)
risk_loci_list = split(risk_loci, cumsum(1:nrow(risk_loci) %in% v))
risk_loci_list = do.call("rbind", lapply(risk_loci_list[2:4], function(y) {return(y[-1,]) } ) )

colnames(risk_loci_list) <- c("SNP", "Chr", "Position", "Closest gene", "Major/minor alleles", "MAF", "Stage1_OR", "Stage1_MetaP", "Stage2_OR","Stage2_MetaP",	"Overall_OR", "Overall_MetaP", "I2_Percent/P")

risk_loci_list$Stage1_MetaP = gsub(" \\× 10","E",risk_loci_list$Stage1_MetaP)
risk_loci_list$Stage1_MetaP = gsub("−","-",risk_loci_list$Stage1_MetaP)

risk_loci_list$Stage2_MetaP = gsub(" \\× 10","E",risk_loci_list$Stage2_MetaP)
risk_loci_list$Stage2_MetaP = gsub("−","-",risk_loci_list$Stage2_MetaP)

risk_loci_list$Overall_MetaP = gsub(" \\× 10","E",risk_loci_list$Overall_MetaP)
risk_loci_list$Overall_MetaP = gsub("−","-",risk_loci_list$Overall_MetaP)

risk_loci_list[,c('Position','MAF','Stage1_MetaP','Stage2_MetaP','Overall_MetaP')] =  sapply(risk_loci_list[,c('Position','MAF','Stage1_MetaP','Stage2_MetaP','Overall_MetaP')], as.numeric )

risk_loci_list$Chr = paste0('chr', risk_loci_list$Chr)
risk_loci_list$chrpos=paste0(risk_loci_list$Chr, ":", risk_loci_list$Position)
risk_loci_list=risk_loci_list[gtools::mixedorder(risk_loci_list$chrpos), ]
risk_loci_list$SNP <- gsub("g","",risk_loci_list$SNP)


## Save results
write.csv(risk_loci_list,file='/dcl01/lieber/ajaffe/Steve/Alz/ad_genetic_risk_loci.csv',row.names=F)
