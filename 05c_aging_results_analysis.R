## Aging analysis
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/aging_controlOnly_mergedStats.rda')

##
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/allRegion_mergedStats_DMP_analysis_dupCor.rda')
load('/dcl01/lieber/ajaffe/Steve/Alz/rdas/merged_DMP_regionSpecific_caseControl_stats.rda')

### Scatterplot of effect sizes?
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_t"])

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_mainEffect_P.Value']),-log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]) )

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_mainEffect_P.Value']),-log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]) )



cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"DLPFC_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"HIPPO_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ERC_Control_Aging_t"])
cor.test(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"CRB_Control_Aging_t"])

###
cor.test(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_P.Value"]) )

cor.test(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]) )

#
png('/dcl01/lieber/ajaffe/Steve/Alz/plots/aging_controlsOnly_vs_AD_Effect.png',res=300,units='in',height=7,width=7)
par(mar=c(5,6,2,2))
plot(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], 
	 Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"ALL_Control_Aging_t"], pch = 21, bg="grey",
	xlab="Cross-region Alzheimer's Effect", ylab="Aging Effect (Controls Only)",
	xlim=c(-12,12), ylim = c(-12,12),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-12,12), y = c(-12,12),col="red")
dev.off()

par(mar=c(5,6,2,2))
plot(-log10(allRegion_mergedStats[,'ALL_subset_interactionEffect_P.Value']), 
	 -log10(Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"INT_Control_Aging_P.Value"]), pch = 21, bg="grey",
	xlab="-log10(Alzheimer's Interaction Pvalue)", ylab="-log10(Aging Interaction Pvalue)",
	xlim=c(0,20), ylim = c(0,20),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(0,20), y = c(0,20),col="red")
dev.off()

table(primary=primaryStats[,'ALL_subset_mainEffect_adj.P.Val']<0.05, 
	 sensitivity=sensitivityStats[rownames(primaryStats),"ALL_subset_mainEffect_adj.P.Val"]<0.05)
