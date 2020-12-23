#####################################
#Make qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/make_qSVs.R)
#####################################

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

## add opioids data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,84, drop=FALSE])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Find genetic race principal components that associate with Group to use in the model below
summary(lm(rse_gene$snpPC1 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC2 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC3 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC8 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC9 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC10 ~ rse_gene$Group))
#These (above) were the components significantly associated with one or both groups (MDD or PTSD)

## model
with(colData(rse_gene), table(Group))

#> with(colData(rse_gene), table(Group))
#Group
#Control     MDD    PTSD 
#    431     432     422 

## also add snpPCs based on association with group
mod = model.matrix(~Group + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10 + opioids_tox, data = colData(rse_gene))

#load cov_rse object from PTSD data
load("rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

#Make sure samples line up across qSVs and gene counts
cov_rse = cov_rse[,colnames(rse_gene)]
#Check
identical(colnames(cov_rse), colnames(rse_gene))

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod) # 19
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

#> getPcaVars(qsvBonf)[1:k]
#[1] 75.500  5.300  3.170  1.760  1.230  1.020  0.816  0.407  0.377  0.352
#[11]  0.326  0.290  0.269  0.242  0.217  0.207  0.189  0.161  0.154

modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_qsvs_sensitivity_opioids.Rdata')
#dir.create("pdf", showWarnings = FALSE)
#pdf(file='../pdf/qsvs_var_explained.pdf', useDingbats = FALSE)
#plot(getPcaVars(qsvBonf)[1:k], pch=20)
#dev.off()



#####################################
#Case-Control BasoAmyg (https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC.R)
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')
library(recount)

#If following from above, don't need to reload objects


###Analysis

#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[c(1, 18:36)] = c("Int", gsub("PC", "qSV", colnames(modQsva)[18:36]))

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

## add opioids data to rse_gene
load("rdas/rse_gene_PTSD_cohort_Feb2019.rda")
colnames(demo)[4] = "BrNum"
rownames(demo) = demo$BrNum
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,84, drop=FALSE])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

#Redefine model using PTSD only
mod = model.matrix(~onlyPTSD + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10 + opioids_tox, data = colData(rse_gene))

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
#> getPcaVars(qsvBonf)[1:k] # 19
#[1] 75.500  5.300  3.170  1.760  1.230  1.020  0.816  0.407  0.377  0.352
#[11]  0.326  0.290  0.269  0.242  0.217  0.207  0.189  0.161  0.154 

modQsva = cbind(mod, qSVs)
save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_onlyPTSD_qsvs_sensitivity_opioids.Rdata')
#pdf(file='../pdf/qsvs_var_explained_onlyPTSD.pdf', useDingbats = FALSE)
#plot(getPcaVars(qsvBonf)[1:k], pch=20)
#dev.off()

#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <- !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]

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
geneStats = cbind(sigGenePTSD[,c(11,13:15)], sigGeneMDD[,c(11,13:15)],sigGeneDx[,c(14:16)],sigGeneonlyPTSD[,c(11,13:15)], PTSDvsMDDPost[,c(11,13:15)])
geneStats = cbind(geneStats, rowData(rse_gene))


## write out qSVA-based stats
geneStats_BasoAmyg = geneStats
save(geneStats_BasoAmyg, file = "rdas/geneStats_DE_qSVA_sensitivity_opioids_BasoAmyg_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_BasoAmyg, file = gzfile("csvs/geneStats_DE_qSVA_sensitivity_opioids_BasoAmyg_threeGroup.csv.gz"))



