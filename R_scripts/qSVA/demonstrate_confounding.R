###

library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(devtools)
library(recount)
library(lattice)

## load counts
load('../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata', verbose = TRUE)

# add MDS (get ethnicity via genotype)
load("../rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# expression filter to remove lowly expressed stuff
#		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

# load qsvBonf, qSVs, mod, modQsva object from PTSD data
load("../rdas/PTSD_qsvs.Rdata", verbose = TRUE)

#filter for MedialAmyg
keepIndex = which(rse_gene$Region == "MedialAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <- !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]

colnames(modQsva)[c(1, 18:36)] = c("Int", gsub("PC", "qSV", colnames(modQsva)[18:36]))

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)

## modeling
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

## get top PTSD
topGenePTSD = topTable(eBGene,coef=3, p.value =0.1,number=nrow(rse_gene))
	
#### which ones are significant
tStatMat = eBGene$t[rownames(topGenePTSD),]
pvalMat = eBGene$p.value[rownames(topGenePTSD),]
rownames(pvalMat) = topGenePTSD$Symbol

###################
## adjusted model #
###################

## modeling
vGeneAdj = voom(dge,mod, plot=TRUE)
fitGeneAdj = lmFit(vGeneAdj)
eBGeneAdj = eBayes(fitGeneAdj)

## get top PTSD
topGenePTSD_Adj = topTable(eBGeneAdj,coef=3, p.value =0.1,number=nrow(rse_gene))[1:20,]
	
#### which ones are significant
tStatMatAdj = eBGeneAdj$t[rownames(topGenePTSD_Adj),]
pvalMatAdj = eBGeneAdj$p.value[rownames(topGenePTSD_Adj),]
rownames(pvalMatAdj) = topGenePTSD_Adj$Symbol

pdf("../pdf/confounding_pvalues_medialAmyg_PTSD.pdf",h=10)
## qSVA on qSVA
pvalMatPlot = pvalMat
pvalMatPlot[pvalMatPlot < 1e-16] = 1e-16
theSeq = seq(0,16,by=0.1)
my.col <- colorRampPalette(c("white","blue"))(length(theSeq))
print(levelplot(-log10(pvalMatPlot[,-1]), aspect = "fill", 
	at = theSeq,pretty=TRUE,xlab="",ylab="",
	main="-log10 P-values: qSVA Model Covariates\namong qSVA Model-significant genes\nPTSD, Medial Amygdala",
	scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
	panel = panel.levelplot.raster, col.regions = my.col))

## Adj on qSVA
pvalMatPlot = eBGene$p.value[rownames(topGenePTSD_Adj),]
rownames(pvalMatPlot) = topGenePTSD_Adj$Symbol
pvalMatPlot[pvalMatPlot < 1e-16] = 1e-16
theSeq = seq(0,16,by=0.1)
my.col <- colorRampPalette(c("white","blue"))(length(theSeq))
print(levelplot(-log10(pvalMatPlot[,-1]), aspect = "fill", 
	at = theSeq,pretty=TRUE,xlab="",ylab="",
	main="-log10 P-values: qSVA Model Covariates\namong Observed Model-significant genes\nPTSD, Medial Amygdala",
	scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
	panel = panel.levelplot.raster, col.regions = my.col))
	
## adj on adj
pvalMatPlot = pvalMatAdj
pvalMatPlot[pvalMatPlot < 1e-16] = 1e-16
theSeq = seq(0,16,by=0.1)
my.col <- colorRampPalette(c("white","blue"))(length(theSeq))
print(levelplot(-log10(pvalMatPlot[,-1]), aspect = "fill", 
	at = theSeq,pretty=TRUE,xlab="",ylab="",
	main="-log10 P-values: Observed Model Covariates\namong Observed Model-significant genes\nPTSD, Medial Amygdala",
	scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
	panel = panel.levelplot.raster, col.regions = my.col))
dev.off()

#####################
## check posthocs ###
#####################

contrast.matrix <- makeContrasts(GroupPTSD, GroupMDD, GroupPTSD-GroupMDD,levels=modQsva)
fit2 = eBayes(contrasts.fit(fitGene, contrast.matrix))

## compare PTSD
ptsdPost = topTable(fit2, coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
ptsdPre = topTable(eBGene,coef=3, p.value =1,sort="none", n =nrow(rse_gene))
all.equal(ptsdPost$t, ptsdPre$t)

## compare MDD
mddPost = topTable(fit2, coef=2,  p.value = 1, sort="none", n = nrow(rse_gene))
mddPre = topTable(eBGene,coef=2, p.value =1,sort="none", n =nrow(rse_gene))
all.equal(mddPost$t, mddPre$t)

## post hoc MDD
mddVsPtsdPost2 = topTable(fit2, coef=3,  p.value = 1, sort="none", n = nrow(rse_gene))

## easy one liner?
posthocContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
mddVsPtsdPost = topTable(eBayes(contrasts.fit(fitGene, posthocContrast)),
	coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
identical(mddVsPtsdPost2, mddVsPtsdPost)