#######
library(jaffelab)
library(edgeR)
library(SummarizedExperiment)
library(RColorBrewer)

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

# #and counts
load("../../../count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")
rse_gene$Group = factor(rse_gene$Group, levels = c("Control", "PTSD", "MDD"))

## read in qSVs from main effects model
load("../../Data/rdas/General/ModelMatrices/Main/EachRegion/PTSD_qsvs.Rdata",
	verbose=TRUE)

## focus on dACC
dIndex =which(rse_gene$Region == "dACC" & 
	rse_gene$PrimaryDx %in% c("Control", "MDD", "PTSD"))
rse_gene_filter = rse_gene[rownames(deg_stats),dIndex]
modQsva_filter = modQsva[colnames(rse_gene_filter),-c(6:8)]

## norm
dge = DGEList(counts = assays(rse_gene_filter)$counts,
	genes = rowData(rse_gene_filter))
dge = calcNormFactors(dge)

#calculate library-size adjustment
vGene = voom(dge,modQsva_filter, plot=FALSE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD vs controls (analysis 1)
outGenePTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")

#### second pass, DEGs within PTSD
bpIndex =which(rse_gene$Region == "dACC" & 
	rse_gene$Group == "PTSD" & 
	rse_gene$PrimaryDx %in% c("Bipolar", "PTSD"))
rse_gene_filter2 = rse_gene[rownames(deg_stats),bpIndex]
modQsva_filter2 = modQsva[colnames(rse_gene_filter2),-c(2:3, 6:8)]
modQsva_filter2 = cbind(modQsva_filter2, 
	PTSD = ifelse(rse_gene_filter2$PrimaryDx == "PTSD", 1, 0))

## norm
dge2 = DGEList(counts = assays(rse_gene_filter2)$counts,
	genes = rowData(rse_gene_filter2))
dge2 = calcNormFactors(dge2)

#calculate library-size adjustment
vGene2 = voom(dge2,modQsva_filter2, plot=FALSE)
fitGene2 = lmFit(vGene2)
eBGene2 = eBayes(fitGene2)

#PTSD vs bipolar (analysis 2)
outGenePTSD2 = topTable(eBGene2,coef=ncol(modQsva_filter2),
	p.value = 1,number=nrow(rse_gene), sort="none")
hist(outGenePTSD2$P.Value)
sum(outGenePTSD2$P.Value < 0.005)
deg_stats[which(outGenePTSD2$adj.P.Val< 0.05) ,]
outGenePTSD2[which(outGenePTSD2$adj.P.Val< 0.05) ,]
min(outGenePTSD2$adj.P.Val) # 0


pdf("figures/figure1/bpd_assessment.pdf")
par(mar=c(5,6,4,2), cex.axis=1.8, cex.lab= 1.8, cex.main=1.8)
plot(outGenePTSD$t, deg_stats$dACC_t_PTSD,
	ylim = c(-6,6),xlim = c(-6,6),
	pch = 21, bg = "grey", main = "dACC: T-stats",
	xlab = "PTSD PrimaryDx (w/o BPD, n=77) vs Control", 
	ylab = "PTSD Group (w/ BPD, n=106) vs Control")
cor(outGenePTSD$t, deg_stats$dACC_t_PTSD) # 0.963
plot(outGenePTSD$logFC, deg_stats$dACC_logFC_PTSD,
	pch = 21, bg = "grey", main = "dACC: log2FC",
	ylim = c(-1.2,1.2),xlim = c(-1.2,1.2),
	xlab = "PTSD PrimaryDx (w/o BPD, n=77) vs Control", 
	ylab = "PTSD Group (w/ BPD, n=106) vs Control")
cor(outGenePTSD$logFC, deg_stats$dACC_logFC_PTSD) # 0.954

plot(outGenePTSD2$t, deg_stats$dACC_t_PTSD,
	ylim = c(-6,6),xlim = c(-6,6),
	pch = 21, bg = "grey", main = "dACC: T-stats",
	xlab = "Primary Dx: PTSD vs BPD (w/in PTSD Group)", 
	ylab = "PTSD Group (w/ BPD, n=106) vs Control")
cor(outGenePTSD2$t, deg_stats$dACC_t_PTSD) # 0.1
dev.off()
