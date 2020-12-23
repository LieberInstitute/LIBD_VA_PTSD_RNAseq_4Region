

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## add manner of death data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,20, drop=FALSE])



## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

## also add snpPCs based on association with group
mod = model.matrix(~Group + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10 + manner_of_death, data = colData(rse_gene))

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
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_qsvs_sensitivity_mannerofdeath.Rdata')


###Analysis

library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')
library(recount)

#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#MDD vs controls (analysis 1)
sigGeneMDD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneMDD) = paste0(colnames(sigGeneMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigGenePTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGenePTSD) = paste0(colnames(sigGenePTSD), "_PTSD")

#Interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")


#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## add manner of death data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,20, drop=FALSE])


## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

#Redefine model using PTSD only
mod = model.matrix(~onlyPTSD + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10 + manner_of_death, data = colData(rse_gene))

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
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]


modQsva = cbind(mod, qSVs)
save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_onlyPTSD_qsvs_sensitivity_mannerofdeath.Rdata')


#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <- !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


###Stuff we want

#Merge qSVA into one big genestats dataframe
geneStats = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx ,sigGeneonlyPTSD, PTSDvsMDDPost)
geneStats = cbind(geneStats, rowData(rse_gene))


## write out qSVA-based stats
geneStats_BasoAmyg = geneStats
save(geneStats_BasoAmyg, file = "rdas/geneStats_DE_qSVA_sensitivity_mannerofdeath_BasoAmyg_threeGroup.rda")
write.csv(geneStats_BasoAmyg, file = gzfile("csvs/geneStats_DE_qSVA_sensitivity_mannerofdeath_BasoAmyg_threeGroup.csv.gz"))

q()

