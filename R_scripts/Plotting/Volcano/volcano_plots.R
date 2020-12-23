#volcano plots
#based on /dcl01/lieber/ajaffe/lab/cst_trap_seq/ip_analysis_genotype/analysis.R
#lines 130 - 138

library(SummarizedExperiment)
library(readxl)
library(jaffelab)
library(edgeR)
library(limma)
library(recount)
library(TeachingDemos) # for shadow text
library(clusterProfiler)
library(org.Mm.eg.db)
library(VariantAnnotation)
library(RColorBrewer)
library(biomaRt)
library(readxl)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")
load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda")
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")

#######
## volanco plots
#MDD
#BLA
geneStats_BasoAmygall$sigColor = as.numeric(geneStats_BasoAmygall$adj.P.Val_MDD < 0.5)+1

pdf("pdf/volcano_plot_BasoAmyg_MDD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_MDD) ~ logFC_MDD, pch = 21, bg = geneStats_BasoAmygall$sigColor,
        data = geneStats_BasoAmygall, xlab = "MDD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "BLA MDD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#PTSD
#BLA
geneStats_BasoAmygall$sigColor = as.numeric(geneStats_BasoAmygall$adj.P.Val_PTSD < 0.5)+1

pdf("pdf/volcano_plot_BasoAmyg_PTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_PTSD) ~ logFC_PTSD, pch = 21, bg = geneStats_BasoAmygall$sigColor,
        data = geneStats_BasoAmygall, xlab = "PTSD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "BLA PTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#onlyPTSD
#BLA
geneStats_BasoAmygall$sigColor = as.numeric(geneStats_BasoAmygall$adj.P.Val_onlyPTSD < 0.5)+1

pdf("pdf/volcano_plot_BasoAmyg_onlyPTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_onlyPTSD) ~ logFC_onlyPTSD, pch = 21, bg = geneStats_BasoAmygall$sigColor,
        data = geneStats_BasoAmygall, xlab = "PTSD vs MDD + Neurotypical log2FC", ylab="-log10(P Value)", main = "BLA onlyPTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()


#MDD
#MeA
geneStats_MedialAmygall$sigColor = as.numeric(geneStats_MedialAmygall$adj.P.Val_MDD < 0.5)+1

pdf("pdf/volcano_plot_MedialAmyg_MDD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_MDD) ~ logFC_MDD, pch = 21, bg = geneStats_MedialAmygall$sigColor,
        data = geneStats_MedialAmygall, xlab = "MDD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "MeA MDD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#PTSD
#MeA
geneStats_MedialAmygall$sigColor = as.numeric(geneStats_MedialAmygall$adj.P.Val_PTSD < 0.5)+1

pdf("pdf/volcano_plot_MedialAmyg_PTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_PTSD) ~ logFC_PTSD, pch = 21, bg = geneStats_MedialAmygall$sigColor,
        data = geneStats_MedialAmygall, xlab = "PTSD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "MeA PTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#onlyPTSD
#MeA
geneStats_MedialAmygall$sigColor = as.numeric(geneStats_MedialAmygall$adj.P.Val_onlyPTSD < 0.5)+1

pdf("pdf/volcano_plot_MedialAmyg_onlyPTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_onlyPTSD) ~ logFC_onlyPTSD, pch = 21, bg = geneStats_MedialAmygall$sigColor,
        data = geneStats_MedialAmygall, xlab = "PTSD vs MDD + Neurotypical log2FC", ylab="-log10(P Value)", main = "MeA onlyPTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()


#MDD
#dACC
geneStats_dACCall$sigColor = as.numeric(geneStats_dACCall$adj.P.Val_MDD < 0.5)+1

pdf("pdf/volcano_plot_dACC_MDD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_MDD) ~ logFC_MDD, pch = 21, bg = geneStats_dACCall$sigColor,
        data = geneStats_dACCall, xlab = "MDD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "dACC MDD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#PTSD
#dACC
geneStats_dACCall$sigColor = as.numeric(geneStats_dACCall$adj.P.Val_PTSD < 0.5)+1

pdf("pdf/volcano_plot_dACC_PTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_PTSD) ~ logFC_PTSD, pch = 21, bg = geneStats_dACCall$sigColor,
        data = geneStats_dACCall, xlab = "PTSD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "dACC PTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#onlyPTSD
#dACC
geneStats_dACCall$sigColor = as.numeric(geneStats_dACCall$adj.P.Val_onlyPTSD < 0.5)+1

pdf("pdf/volcano_plot_dACC_onlyPTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_onlyPTSD) ~ logFC_onlyPTSD, pch = 21, bg = geneStats_dACCall$sigColor,
        data = geneStats_dACCall, xlab = "PTSD vs MDD + Neurotypical log2FC", ylab="-log10(P Value)", main = "dACC onlyPTSD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()


#MDD
#DLPFC
geneStats_DLPFCall$sigColor = as.numeric(geneStats_DLPFCall$adj.P.Val_MDD < 0.5)+1

pdf("pdf/volcano_plot_DLPFC_MDD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_MDD) ~ logFC_MDD, pch = 21, bg = geneStats_DLPFCall$sigColor,
        data = geneStats_DLPFCall, xlab = "MDD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "DLPFC MDD DEGs")
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#PTSD
#DLPFC
geneStats_DLPFCall <- geneStats_DLPFCall[order(geneStats_DLPFCall$P.Value_PTSD),]

geneStats_DLPFCall$sigColor = as.numeric(geneStats_DLPFCall$adj.P.Val_PTSD < 0.1)+1

pdf("pdf/volcano_plot_DLPFC_PTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_PTSD) ~ logFC_PTSD, pch = 21, bg = geneStats_DLPFCall$sigColor,
        data = geneStats_DLPFCall, xlab = "PTSD vs Neurotypical log2FC", ylab="-log10(P Value)", main = "DLPFC PTSD DEGs", ylim = c(0,6))
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()

#onlyPTSD
#DLPFC
geneStats_DLPFCall <- geneStats_DLPFCall[order(geneStats_DLPFCall$P.Value_onlyPTSD),]
geneStats_DLPFCall$sigColor = as.numeric(geneStats_DLPFCall$adj.P.Val_onlyPTSD < 0.5)+1

pdf("pdf/volcano_plot_DLPFC_onlyPTSD.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Set1"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value_onlyPTSD) ~ logFC_onlyPTSD, pch = 21, bg = geneStats_DLPFCall$sigColor,
        data = geneStats_DLPFCall, xlab = "PTSD vs MDD + Neurotypical log2FC", ylab="-log10(P Value)", main = "DLPFC onlyPTSD DEGs")
text(-log10(geneStats_DLPFCall$P.Value_onlyPTSD)[1:6] ~ geneStats_DLPFCall$logFC_onlyPTSD[1:6], labels=geneStats_DLPFCall$Symbol[1:6], data=geneStats_DLPFCall[1:6], cex=0.9, font=2)
abline(h=-log10(0.05), lty=2,lwd=2)
dev.off()






#scatter
DLPFC <- summary(lm(geneStats_DLPFCall$t_onlyPTSD~geneStats_DLPFCall$t_PTSD))

pdf("pdf/tstat_Scatterplot_DLPFC_PTSD.pdf",useDingbats=FALSE)
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(t_onlyPTSD ~ t_PTSD, pch = 21, bg = "#377eb8",
	 data = geneStats_DLPFCall, xlab = "t-statistics PTSD", ylab="t-statistics onlyPTSD", main = "DLPFC t-Statistics Comparison")
legend("topleft", bty="n", legend=paste("R2 = ", format(DLPFC$adj.r.squared, digits=3)))
dev.off()




#scatter
pdf("pdf/FC_Scatterplot_DLPFC_PTSD.pdf",useDingbats=FALSE)
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(logFC_onlyPTSD ~ logFC_PTSD, pch = 21, bg = "#377eb8",
	 data = geneStats_DLPFCall, xlab = "log2FC PTSD", ylab="log2FC onlyPTSD", main = "DLPFC log2FC Comparison")
dev.off()



#scatter
pdf("pdf/FC_Scatterplot_DLPFC_MDDonlyPTSD.pdf",useDingbats=FALSE)
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(logFC_onlyPTSD ~ logFC_MDD, pch = 21, bg = "#377eb8",
	 data = geneStats_DLPFCall, xlab = "log2FC MDD", ylab="log2FC onlyPTSD", main = "DLPFC log2FC Comparison")
dev.off()

#scatter
DLPFC2 <- summary(lm(geneStats_DLPFCall$t_onlyPTSD~geneStats_DLPFCall$t_MDD))

pdf("pdf/tstat_Scatterplot_DLPFC_MDDonlyPTSD.pdf",useDingbats=FALSE)
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(t_onlyPTSD ~ t_MDD, pch = 21, bg = "#377eb8",
	 data = geneStats_DLPFCall, xlab = "t-statistics MDD", ylab="t-statistics onlyPTSD", main = "DLPFC t-Statistics Comparison")
legend("topleft", bty="n", legend=paste("R2 = ", format(DLPFC2$adj.r.squared, digits=3)))
dev.off()



geneStats_BasoAmygall <- geneStats_BasoAmygall[order(geneStats_BasoAmygall$P.Value_MDD),]
