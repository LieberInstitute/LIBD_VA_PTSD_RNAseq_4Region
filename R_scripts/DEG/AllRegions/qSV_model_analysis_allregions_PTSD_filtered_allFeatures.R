
############################################################################################33
#Load up

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])
#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])
#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])
#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])


## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]
#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]
#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]
#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]

#Remake model for across region analyses
mod = model.matrix(~Group*Region + AgeDeath + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))
#rename colnames to be syntactically valid
colnames(mod) <- make.names(colnames(mod))
colnames(mod)[1] = "Int"

#load cov_rse object from PTSD data
load("rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

#Make sure samples line up across qSVs and gene counts
cov_rse = cov_rse[,colnames(rse_gene)]
#Check
identical(colnames(cov_rse), colnames(rse_gene))

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
k
qSV_mat = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

modQsva = cbind(mod, qSV_mat)

save(qsvBonf, qSV_mat, modQsva, mod, file = 'rdas/PTSD_qsvs_Regionintxn.Rdata')

##### GENE ######

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=FALSE)

## do duplicate correlation
gene_dupCorr = duplicateCorrelation(vGene$E, modQsva, block=colData(rse_gene)$BrNum)
save(gene_dupCorr, file = "rdas/geneLevel_duplicateCorrelation_limma_forDE_PTSD.rda")

# and then fit
fitGene = lmFit(vGene, modQsva,
        correlation=gene_dupCorr$consensus.correlation,
        block=colData(rse_gene)$BrNum)
eBGene = eBayes(fitGene)

#MDD vs controls (analysis 1)
sigGeneMDD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneMDD) = paste0(colnames(sigGeneMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigGenePTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGenePTSD) = paste0(colnames(sigGenePTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostGene = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPostGene) = paste0(colnames(PTSDvsMDDPostGene), "_PTSDvsMDD")

#Region interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=c(21,23,25),
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_rintxnMDD")

#Region interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=c(22,24,26),
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigGeneRegion = topTable(eBGene,coef=c(4,5,6),
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneRegion) = paste0(colnames(sigGeneRegion), "_Region")

												  
##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
        genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)

## do duplicate correlation
exon_dupCorr = duplicateCorrelation(vExon$E, modQsva, block=colData(rse_exon)$BrNum)
save(exon_dupCorr, file = "rdas/exonLevel_duplicateCorrelation_limma_forDE_PTSD.rda")

# and then fit
fitExon = lmFit(vExon, modQsva,
        correlation=exon_dupCorr$consensus.correlation,
        block=colData(rse_exon)$BrNum)
eBExon = eBayes(fitExon)

#MDD vs controls (analysis 1)
sigExonMDD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonMDD) = paste0(colnames(sigExonMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigExonPTSD = topTable(eBExon,coef=3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonPTSD) = paste0(colnames(sigExonPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigExonDx = topTable(eBExon,coef=2:3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonDx) = paste0(colnames(sigExonDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostExon = topTable(eBayes(contrasts.fit(fitExon, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_exon))
colnames(PTSDvsMDDPostExon) = paste0(colnames(PTSDvsMDDPostExon), "_PTSDvsMDD")

#Region interaction MDD
outExon_interactionEffect_MDD = topTable(eBExon,coef=c(21,23,25),
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_MDD) = paste0(colnames(outExon_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outExon_interactionEffect_PTSD = topTable(eBExon,coef=c(22,24,26),
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_PTSD) = paste0(colnames(outExon_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigExonRegion = topTable(eBExon,coef=c(4,5,6),
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonRegion) = paste0(colnames(sigExonRegion), "_Region")


##### JUNCTION ######

dje = DGEList(counts = assays(rse_jx)$counts,
        genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)

## do duplicate correlation
jx_dupCorr = duplicateCorrelation(vJxn$E, modQsva, block=colData(rse_jx)$BrNum)
save(jx_dupCorr, file = "rdas/jxLevel_duplicateCorrelation_limma_forDE_PTSD.rda")

# and then fit
fitJxn = lmFit(vJxn, modQsva,
        correlation=jx_dupCorr$consensus.correlation,
        block=colData(rse_jx)$BrNum)
eBJxn = eBayes(fitJxn)

#MDD vs controls (analysis 1)
sigJxnMDD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnMDD) = paste0(colnames(sigJxnMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigJxnPTSD = topTable(eBJxn,coef=3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnPTSD) = paste0(colnames(sigJxnPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigJxnDx = topTable(eBJxn,coef=2:3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnDx) = paste0(colnames(sigJxnDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostJxn = topTable(eBayes(contrasts.fit(fitJxn, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_jx))
colnames(PTSDvsMDDPostJxn) = paste0(colnames(PTSDvsMDDPostJxn), "_PTSDvsMDD")

#Region interaction MDD
outJxn_interactionEffect_MDD = topTable(eBJxn,coef=c(21,23,25),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_MDD) = paste0(colnames(outJxn_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outJxn_interactionEffect_PTSD = topTable(eBJxn,coef=c(22,24,26),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_PTSD) = paste0(colnames(outJxn_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigJxnRegion = topTable(eBJxn,coef=c(4,5,6),
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnRegion) = paste0(colnames(sigJxnRegion), "_Region")


##### TRANSCRIPT ######

txExprs = log2(assays(rse_tx)$tpm+ 1)

## do duplicate correlation
tx_dupCorr = duplicateCorrelation(txExprs, modQsva, block=colData(rse_tx)$BrNum)
save(tx_dupCorr, file = "rdas/txLevel_duplicateCorrelation_limma_forDE_PTSD.rda")

# and then fit
fitTx = lmFit(txExprs, modQsva,
        correlation=tx_dupCorr$consensus.correlation,
        block=colData(rse_tx)$BrNum)
eBTx = eBayes(fitTx)

#MDD vs controls (analysis 1)
sigTxMDD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxMDD) = paste0(colnames(sigTxMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigTxPTSD = topTable(eBTx,coef=3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxPTSD) = paste0(colnames(sigTxPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigTxDx = topTable(eBTx,coef=2:3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxDx) = paste0(colnames(sigTxDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostTx = topTable(eBayes(contrasts.fit(fitTx, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, genelist = rowRanges(rse_tx), sort="none", n = nrow(rse_tx))
colnames(PTSDvsMDDPostTx) = paste0(colnames(PTSDvsMDDPostTx), "_PTSDvsMDD")

#Region interaction MDD
outTx_interactionEffect_MDD = topTable(eBTx,coef=c(21,23,25),
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_MDD) = paste0(colnames(outTx_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outTx_interactionEffect_PTSD = topTable(eBTx,coef=c(22,24,26),
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_PTSD) = paste0(colnames(outTx_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigTxRegion = topTable(eBTx,coef=c(4,5,6),
	p.value = 1,number=nrow(rse_tx), sort="none")
colnames(sigTxRegion) = paste0(colnames(sigTxRegion), "_Region")

#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#reload rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])
#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])
#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])
#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]
#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]
#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]
#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 
rse_exon$onlyPTSD = ifelse(rse_exon$Group == "PTSD" ,1, 0) 
rse_jx$onlyPTSD = ifelse(rse_jx$Group == "PTSD" ,1, 0) 
rse_tx$onlyPTSD = ifelse(rse_tx$Group == "PTSD" ,1, 0) 

#Redefine model using PTSD only
mod = model.matrix(~onlyPTSD*Region + AgeDeath + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))
#rename colnames to be syntactically valid
colnames(mod) <- make.names(colnames(mod))
colnames(mod)[1] = "Int"

#load cov_rse object from PTSD data
load("rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

#Make sure samples line up across qSVs and gene counts
cov_rse = cov_rse[,colnames(rse_gene)]
#Check
identical(colnames(cov_rse), colnames(rse_gene))

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod) 
k
qSV_mat = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

modQsva = cbind(mod, qSV_mat)
save(qsvBonf, qSV_mat, modQsva, mod, file = 'rdas/PTSD_qsvs_Regionintxn_onlyPTSD.Rdata')

##### GENE ######

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)

## do duplicate correlation
gene_dupCorr = duplicateCorrelation(vGene$E, modQsva, block=colData(rse_gene)$BrNum)
save(gene_dupCorr, file = "rdas/geneLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda")

# and then fit
fitGene = lmFit(vGene, modQsva,
        correlation=gene_dupCorr$consensus.correlation,
        block=colData(rse_gene)$BrNum)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=c(20:22),
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")

#Region effect
sigGeneRegion_onlyPTSD = topTable(eBGene,coef=c(3,4,5),
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneRegion_onlyPTSD) = paste0(colnames(sigGeneRegion_onlyPTSD), "_Region_onlyPTSD")

##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
        genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)

## do duplicate correlation
exon_dupCorr = duplicateCorrelation(vExon$E, modQsva, block=colData(rse_exon)$BrNum)
save(exon_dupCorr, file = "rdas/exonLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda")

# and then fit
fitExon = lmFit(vExon, modQsva,
        correlation=exon_dupCorr$consensus.correlation,
        block=colData(rse_exon)$BrNum)
eBExon = eBayes(fitExon)

#PTSD only
sigExononlyPTSD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExononlyPTSD) = paste0(colnames(sigExononlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outExon_interactionEffect_onlyPTSD = topTable(eBExon,coef=c(20:22),
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_onlyPTSD) = paste0(colnames(outExon_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")

#Region effect
sigExonRegion_onlyPTSD = topTable(eBExon,coef=c(3,4,5),
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonRegion_onlyPTSD) = paste0(colnames(sigExonRegion_onlyPTSD), "_Region_onlyPTSD")


##### JUNCTION ######

dje = DGEList(counts = assays(rse_jx)$counts,
        genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)

## do duplicate correlation
jx_dupCorr = duplicateCorrelation(vJxn$E, modQsva, block=colData(rse_jx)$BrNum)
save(jx_dupCorr, file = "rdas/jxLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda")

# and then fit
fitJxn = lmFit(vJxn, modQsva,
        correlation=jx_dupCorr$consensus.correlation,
        block=colData(rse_jx)$BrNum)
eBJxn = eBayes(fitJxn)

#PTSD only
sigJxnonlyPTSD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnonlyPTSD) = paste0(colnames(sigJxnonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outJxn_interactionEffect_onlyPTSD = topTable(eBJxn,coef=c(20:22),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_onlyPTSD) = paste0(colnames(outJxn_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")
													  
#Region effect
sigJxnRegion_onlyPTSD = topTable(eBJxn,coef=c(3,4,5),
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnRegion_onlyPTSD) = paste0(colnames(sigJxnRegion_onlyPTSD), "_Region_onlyPTSD")
													  
##### TRANSCRIPT ######

txExprs = log2(assays(rse_tx)$tpm+ 1)

## do duplicate correlation
tx_dupCorr = duplicateCorrelation(txExprs, modQsva, block=colData(rse_tx)$BrNum)
save(tx_dupCorr, file = "rdas/txLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda")

# and then fit
fitTx = lmFit(txExprs, modQsva,
        correlation=tx_dupCorr$consensus.correlation,
        block=colData(rse_tx)$BrNum)
eBTx = eBayes(fitTx)

#PTSD only
sigTxonlyPTSD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxonlyPTSD) = paste0(colnames(sigTxonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outTx_interactionEffect_onlyPTSD = topTable(eBTx,coef=c(20:22),
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_onlyPTSD) = paste0(colnames(outTx_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")
								
#Region effect
sigTxRegion_onlyPTSD = topTable(eBTx,coef=c(3,4,5),
	p.value = 1,number=nrow(rse_tx), sort="none")
colnames(sigTxRegion_onlyPTSD) = paste0(colnames(sigTxRegion_onlyPTSD), "_Region_onlyPTSD")
												
###Stuff we want

###### GENE ########
#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPostGene, sigGeneRegion, outGene_interactionEffect_MDD, 
	outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD, sigGeneRegion_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda")

###### EXON ########
#All columns
exonStatsall = cbind(sigExonPTSD, sigExonMDD, sigExonDx, sigExononlyPTSD, PTSDvsMDDPostExon, sigExonRegion, outExon_interactionEffect_MDD, 
	outExon_interactionEffect_PTSD, outExon_interactionEffect_onlyPTSD,sigExonRegion_onlyPTSD)
exonStatsall = cbind(exonStatsall, rowData(rse_exon))
## write out qSVA-based stats
save(exonStatsall, file = "rdas/all_regions/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda")

###### JUNCTION ########
#All columns
jxnStatsall = cbind(sigJxnPTSD, sigJxnMDD, sigJxnDx, sigJxnonlyPTSD, PTSDvsMDDPostJxn, sigJxnRegion, outJxn_interactionEffect_MDD, 
	outJxn_interactionEffect_PTSD, outJxn_interactionEffect_onlyPTSD,sigJxnRegion_onlyPTSD)
jxnStatsall = cbind(jxnStatsall, rowData(rse_jx))
## write out qSVA-based stats
save(jxnStatsall, file = "rdas/all_regions/jxStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda")

####### TRANSCRIPT ########
#All columns
txStatsall = cbind(sigTxPTSD, sigTxMDD, sigTxDx, sigTxonlyPTSD, PTSDvsMDDPostTx, sigTxRegion, outTx_interactionEffect_MDD, 
	outTx_interactionEffect_PTSD, outTx_interactionEffect_onlyPTSD,sigTxRegion_onlyPTSD)
txStatsall = cbind(txStatsall, rowData(rse_tx))
## write out qSVA-based stats
save(txStatsall, file = "rdas/all_regions/txStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda")



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()
