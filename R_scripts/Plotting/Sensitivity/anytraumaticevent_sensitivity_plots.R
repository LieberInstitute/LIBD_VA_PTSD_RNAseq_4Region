###################Get R2 values######################

#load everything
load('rdas/dACC/geneStats_DE_qSVA_sensitivity_any_traumatic_event_dACC_threeGroup.rda')
dACC_any_traumatic_event <- geneStats_dACC
load('rdas/dACC/geneStats_DE_qSVA_dACC_threeGroup.rda')

load('rdas/DLPFC/geneStats_DE_qSVA_sensitivity_any_traumatic_event_DLPFC_threeGroup.rda')
DLPFC_any_traumatic_event <- geneStats_DLPFC
load('rdas/DLPFC/geneStats_DE_qSVA_DLPFC_threeGroup.rda')

load('rdas/BasoAmyg/geneStats_DE_qSVA_sensitivity_any_traumatic_event_BasoAmyg_threeGroup.rda')
BasoAmyg_any_traumatic_event <- geneStats_BasoAmyg
load('rdas/BasoAmyg/geneStats_DE_qSVA_BasoAmyg_threeGroup.rda')

load('rdas/MedialAmyg/geneStats_DE_qSVA_sensitivity_any_traumatic_event_MedialAmyg_threeGroup.rda')
MedialAmyg_any_traumatic_event <- geneStats_MedialAmyg
load('rdas/MedialAmyg/geneStats_DE_qSVA_MedialAmyg_threeGroup.rda')
geneStats_MedialAmyg <- geneStats_medialamyg

#################################onlyPTSD###################################

dACC_onlyPTSD <- summary(lm(dACC_any_traumatic_event$logFC_onlyPTSD~geneStats_dACC$logFC_onlyPTSD))
DLPFC_onlyPTSD <- summary(lm(DLPFC_any_traumatic_event$logFC_onlyPTSD~geneStats_DLPFC$logFC_onlyPTSD))
BasoAmyg_onlyPTSD <- summary(lm(BasoAmyg_any_traumatic_event$logFC_onlyPTSD~geneStats_BasoAmyg$logFC_onlyPTSD))
MedialAmyg_onlyPTSD <- summary(lm(MedialAmyg_any_traumatic_event$logFC_onlyPTSD~geneStats_MedialAmyg$logFC_onlyPTSD))


#################################PTSD###################################


dACC_PTSD <- summary(lm(dACC_any_traumatic_event$logFC_PTSD~geneStats_dACC$logFC_PTSD))
DLPFC_PTSD <- summary(lm(DLPFC_any_traumatic_event$logFC_PTSD~geneStats_DLPFC$logFC_PTSD))
BasoAmyg_PTSD <- summary(lm(BasoAmyg_any_traumatic_event$logFC_PTSD~geneStats_BasoAmyg$logFC_PTSD))
MedialAmyg_PTSD <- summary(lm(MedialAmyg_any_traumatic_event$logFC_PTSD~geneStats_MedialAmyg$logFC_PTSD))


#################################MDD###################################


dACC_MDD <- summary(lm(dACC_any_traumatic_event$logFC_MDD~geneStats_dACC$logFC_MDD))
DLPFC_MDD <- summary(lm(DLPFC_any_traumatic_event$logFC_MDD~geneStats_DLPFC$logFC_MDD))
BasoAmyg_MDD <- summary(lm(BasoAmyg_any_traumatic_event$logFC_MDD~geneStats_BasoAmyg$logFC_MDD))
MedialAmyg_MDD <- summary(lm(MedialAmyg_any_traumatic_event$logFC_MDD~geneStats_MedialAmyg$logFC_MDD))

#################################PTSDvsMDD###############################


dACC_PTSDvsMDD <- summary(lm(dACC_any_traumatic_event$logFC_PTSDvsMDD~geneStats_dACC$logFC_PTSDvsMDD))
DLPFC_PTSDvsMDD <- summary(lm(DLPFC_any_traumatic_event$logFC_PTSDvsMDD~geneStats_DLPFC$logFC_PTSDvsMDD))
BasoAmyg_PTSDvsMDD <- summary(lm(BasoAmyg_any_traumatic_event$logFC_PTSDvsMDD~geneStats_BasoAmyg$logFC_PTSDvsMDD))
MedialAmyg_PTSDvsMDD <- summary(lm(MedialAmyg_any_traumatic_event$logFC_PTSDvsMDD~geneStats_MedialAmyg$logFC_PTSDvsMDD))

#################################ANOVA###################################


dACC_ANOVA <- summary(lm(dACC_any_traumatic_event$F_ANOVA~geneStats_dACC$F_ANOVA))
DLPFC_ANOVA <- summary(lm(DLPFC_any_traumatic_event$F_ANOVA~geneStats_DLPFC$F_ANOVA))
BasoAmyg_ANOVA <- summary(lm(BasoAmyg_any_traumatic_event$F_ANOVA~geneStats_BasoAmyg$F_ANOVA))
MedialAmyg_ANOVA <- summary(lm(MedialAmyg_any_traumatic_event$F_ANOVA~geneStats_MedialAmyg$F_ANOVA))




###############################
#########Plots#################
###############################


#################################onlyPTSD###################################

pdf('pdf/Sensitivity_Analyses/any_traumatic_event_sensitivity_plots.pdf',useDingbats = FALSE)
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(dACC_any_traumatic_event$logFC_onlyPTSD, geneStats_dACC$logFC_onlyPTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "dACC")
abline(h=mean(geneStats_dACC$logFC_onlyPTSD), v=mean(dACC_any_traumatic_event$logFC_onlyPTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(dACC_onlyPTSD$adj.r.squared, digits=3)))
plot(DLPFC_any_traumatic_event$logFC_onlyPTSD, geneStats_DLPFC$logFC_onlyPTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "DLPFC")
abline(h=mean(geneStats_DLPFC$logFC_onlyPTSD), v=mean(DLPFC_any_traumatic_event$logFC_onlyPTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(DLPFC_onlyPTSD$adj.r.squared, digits=3)))
plot(BasoAmyg_any_traumatic_event$logFC_onlyPTSD, geneStats_BasoAmyg$logFC_onlyPTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "BasoAmyg")
abline(h=mean(geneStats_BasoAmyg$logFC_onlyPTSD), v=mean(BasoAmyg_any_traumatic_event$logFC_onlyPTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(BasoAmyg_onlyPTSD$adj.r.squared, digits=3)))
plot(MedialAmyg_any_traumatic_event$logFC_onlyPTSD, geneStats_MedialAmyg$logFC_onlyPTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "MedialAmyg")
abline(h=mean(geneStats_MedialAmyg$logFC_onlyPTSD), v=mean(MedialAmyg_any_traumatic_event$logFC_onlyPTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(MedialAmyg_onlyPTSD$adj.r.squared, digits=3)))
mtext("Gene Expression Comparison - Any Traumatic Event Sensitivity Analysis onlyPTSD", outer = TRUE, cex = 1.0)

#################################PTSD###################################


par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(dACC_any_traumatic_event$logFC_PTSD, geneStats_dACC$logFC_PTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "dACC")
abline(h=mean(geneStats_dACC$logFC_PTSD), v=mean(dACC_any_traumatic_event$logFC_PTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(dACC_PTSD$adj.r.squared, digits=3)))
plot(DLPFC_any_traumatic_event$logFC_PTSD, geneStats_DLPFC$logFC_PTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "DLPFC")
abline(h=mean(geneStats_DLPFC$logFC_PTSD), v=mean(DLPFC_any_traumatic_event$logFC_PTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(DLPFC_PTSD$adj.r.squared, digits=3)))
plot(BasoAmyg_any_traumatic_event$logFC_PTSD, geneStats_BasoAmyg$logFC_PTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "BasoAmyg")
abline(h=mean(geneStats_BasoAmyg$logFC_PTSD), v=mean(BasoAmyg_any_traumatic_event$logFC_PTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(BasoAmyg_PTSD$adj.r.squared, digits=3)))
plot(MedialAmyg_any_traumatic_event$logFC_PTSD, geneStats_MedialAmyg$logFC_PTSD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "MedialAmyg")
abline(h=mean(geneStats_MedialAmyg$logFC_PTSD), v=mean(MedialAmyg_any_traumatic_event$logFC_PTSD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(MedialAmyg_PTSD$adj.r.squared, digits=3)))
mtext("Gene Expression Comparison - Any Traumatic Event Sensitivity Analysis PTSD", outer = TRUE, cex = 1.0)

#################################MDD###################################


par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(dACC_any_traumatic_event$logFC_MDD, geneStats_dACC$logFC_MDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "dACC")
abline(h=mean(geneStats_dACC$logFC_MDD), v=mean(dACC_any_traumatic_event$logFC_MDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(dACC_MDD$adj.r.squared, digits=3)))
plot(DLPFC_any_traumatic_event$logFC_MDD, geneStats_DLPFC$logFC_MDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "DLPFC")
abline(h=mean(geneStats_DLPFC$logFC_MDD), v=mean(DLPFC_any_traumatic_event$logFC_MDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(DLPFC_MDD$adj.r.squared, digits=3)))
plot(BasoAmyg_any_traumatic_event$logFC_MDD, geneStats_BasoAmyg$logFC_MDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "BasoAmyg")
abline(h=mean(geneStats_BasoAmyg$logFC_MDD), v=mean(BasoAmyg_any_traumatic_event$logFC_MDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(BasoAmyg_MDD$adj.r.squared, digits=3)))
plot(MedialAmyg_any_traumatic_event$logFC_MDD, geneStats_MedialAmyg$logFC_MDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "MedialAmyg")
abline(h=mean(geneStats_MedialAmyg$logFC_MDD), v=mean(MedialAmyg_any_traumatic_event$logFC_MDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(MedialAmyg_MDD$adj.r.squared, digits=3)))
mtext("Gene Expression Comparison - Any Traumatic Event Sensitivity Analysis MDD", outer = TRUE, cex = 1.0)

#################################PTSDvsMDD###################################


par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(dACC_any_traumatic_event$logFC_PTSDvsMDD, geneStats_dACC$logFC_PTSDvsMDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "dACC")
abline(h=mean(geneStats_dACC$logFC_PTSDvsMDD), v=mean(dACC_any_traumatic_event$logFC_PTSDvsMDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(dACC_PTSDvsMDD$adj.r.squared, digits=3)))
plot(DLPFC_any_traumatic_event$logFC_PTSDvsMDD, geneStats_DLPFC$logFC_PTSDvsMDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "DLPFC")
abline(h=mean(geneStats_DLPFC$logFC_PTSDvsMDD), v=mean(DLPFC_any_traumatic_event$logFC_PTSDvsMDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(DLPFC_PTSDvsMDD$adj.r.squared, digits=3)))
plot(BasoAmyg_any_traumatic_event$logFC_PTSDvsMDD, geneStats_BasoAmyg$logFC_PTSDvsMDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "BasoAmyg")
abline(h=mean(geneStats_BasoAmyg$logFC_PTSDvsMDD), v=mean(BasoAmyg_any_traumatic_event$logFC_PTSDvsMDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(BasoAmyg_PTSDvsMDD$adj.r.squared, digits=3)))
plot(MedialAmyg_any_traumatic_event$logFC_PTSDvsMDD, geneStats_MedialAmyg$logFC_PTSDvsMDD, xlab = "Including Any Traumatic Event log2FC", ylab="Original log2FC", main = "MedialAmyg")
abline(h=mean(geneStats_MedialAmyg$logFC_PTSDvsMDD), v=mean(MedialAmyg_any_traumatic_event$logFC_PTSDvsMDD), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(MedialAmyg_PTSDvsMDD$adj.r.squared, digits=3)))
mtext("Gene Expression Comparison - Any Traumatic Event Sensitivity Analysis PTSDvsMDD", outer = TRUE, cex = 1.0)

#################################ANOVA###################################


par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(dACC_any_traumatic_event$F_ANOVA, geneStats_dACC$F_ANOVA, xlab = "Including Any Traumatic Event F-Statistics", ylab="Original F-Statistics", main = "dACC")
abline(h=mean(geneStats_dACC$F_ANOVA), v=mean(dACC_any_traumatic_event$F_ANOVA), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(dACC_ANOVA$adj.r.squared, digits=3)))
plot(DLPFC_any_traumatic_event$F_ANOVA, geneStats_DLPFC$F_ANOVA, xlab = "Including Any Traumatic Event F-Statistics", ylab="Original F-Statistics", main = "DLPFC")
abline(h=mean(geneStats_DLPFC$F_ANOVA), v=mean(DLPFC_any_traumatic_event$F_ANOVA), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(DLPFC_ANOVA$adj.r.squared, digits=3)))
plot(BasoAmyg_any_traumatic_event$F_ANOVA, geneStats_BasoAmyg$F_ANOVA, xlab = "Including Any Traumatic Event F-Statistics", ylab="Original F-Statistics", main = "BasoAmyg")
abline(h=mean(geneStats_BasoAmyg$F_ANOVA), v=mean(BasoAmyg_any_traumatic_event$F_ANOVA), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(BasoAmyg_ANOVA$adj.r.squared, digits=3)))
plot(MedialAmyg_any_traumatic_event$F_ANOVA, geneStats_MedialAmyg$F_ANOVA, xlab = "Including Any Traumatic Event F-Statistics", ylab="Original F-Statistics", main = "MedialAmyg")
abline(h=mean(geneStats_MedialAmyg$F_ANOVA), v=mean(MedialAmyg_any_traumatic_event$F_ANOVA), col="purple")
legend("bottomright", bty="n", legend=paste("R2 = ", format(MedialAmyg_ANOVA$adj.r.squared, digits=3)))
mtext("Gene Expression Comparison - Any Traumatic Event Sensitivity Analysis ANOVA", outer = TRUE, cex = 1.0)
dev.off()
