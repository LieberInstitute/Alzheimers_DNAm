### Check old stats (n380) versus new (n377)
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/caseControl_DMC_allRegion.rda')
oldStats=allStats
oldStats_sig = oldStats$Name[oldStats$Primary_subset_mainEffect_adj.P.Val<0.05]

max(oldStats[oldStats$Primary_subset_mainEffect_adj.P.Val<0.05,'Primary_subset_mainEffect_P.Value'])

## new stats
load('/dcl01/lieber/ajaffe/Steve/Alz/Paper/rdas/caseControl_DMC_allRegion.rda')
newStats=allStats
newStats_sig = newStats$Name[newStats$Primary_subset_mainEffect_adj.P.Val<0.05]
max(newStats[newStats$Primary_subset_mainEffect_adj.P.Val<0.05,'Primary_subset_mainEffect_P.Value'])


table(oldSig=oldStats$Primary_subset_mainEffect_adj.P.Val<0.05, newSig=oldStats$Name %in% newStats_sig )
##
cor.test(oldStats$Primary_subset_mainEffect_logFC,newStats[oldStats$Name,'Primary_subset_mainEffect_logFC'])

plot(oldStats$Primary_subset_mainEffect_logFC,newStats[oldStats$Name,'Primary_subset_mainEffect_logFC'])

plot(oldStats[oldStats_sig,'Primary_subset_mainEffect_logFC'],
newStats[oldStats_sig,'Primary_subset_mainEffect_logFC'])

plot(-log10(oldStats[oldStats_sig,'Primary_subset_mainEffect_P.Value']),
-log10(newStats[oldStats_sig,'Primary_subset_mainEffect_P.Value']) )


plot(oldStats[newStats_sig,'Primary_subset_mainEffect_logFC'],
newStats[newStats_sig,'Primary_subset_mainEffect_logFC'])

plot(-log10(oldStats[newStats_sig,'Primary_subset_mainEffect_P.Value']),
-log10(newStats[newStats_sig,'Primary_subset_mainEffect_P.Value']) )

plot(-log10(oldStats[oldStats_sig[!oldStats_sig%in%newStats_sig],'Primary_subset_mainEffect_P.Value']),
-log10(newStats[oldStats_sig[!oldStats_sig%in%newStats_sig],'Primary_subset_mainEffect_P.Value']) )
abline(0,1)


plot(oldStats[oldStats_sig[!oldStats_sig%in%newStats_sig],'Primary_subset_mainEffect_logFC'],
newStats[oldStats_sig[!oldStats_sig%in%newStats_sig],'Primary_subset_mainEffect_logFC'])
abline(0,1)
