
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

## add lifetime antidepressants data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
#change encoding of lifetime_antidepress to be 0 or 1 instead of true or false
#Note, changes the 9 "NA" samples to 0
demo$lifetime_antidepress <- ifelse(demo$lifetime_antidepress == "TRUE",1,0)
demo$lifetime_antidepress[is.na(demo$lifetime_antidepress)] <- 0
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,79, drop=FALSE])



## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

load('rdas/PTSD_qsvs_sensitivity_lifetimeantidepressants.Rdata')


###Analysis
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')
library(recount)

#filter for MedialAmyg 
keepIndex = which(rse_gene$Region == "MedialAmyg")
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

## add lifetime antidepressants data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
#change encoding of lifetime_antidepress to be 0 or 1 instead of true or false
#Note, changes the 9 "NA" samples to 0
demo$lifetime_antidepress <- ifelse(demo$lifetime_antidepress == "TRUE",1,0)
demo$lifetime_antidepress[is.na(demo$lifetime_antidepress)] <- 0
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,79, drop=FALSE])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_sensitivity_lifetimeantidepressants.Rdata')


#filter for MedialAmyg 
keepIndex = which(rse_gene$Region == "MedialAmyg")
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
geneStats_MedialAmyg = geneStats
save(geneStats_MedialAmyg, file = "rdas/geneStats_DE_qSVA_sensitivity_lifetimeantidepressants_MedialAmyg_threeGroup.rda")
write.csv(geneStats_MedialAmyg, file = gzfile("csvs/geneStats_DE_qSVA_sensitivity_lifetimeantidepressants_MedialAmyg_threeGroup.csv.gz"))


q()

