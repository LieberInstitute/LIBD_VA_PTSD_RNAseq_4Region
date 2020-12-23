#####Stratify by dx analysis#######


############# Controls/MDD ##############

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')




#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

#Redefine model using PTSD only
mod = model.matrix(~Sex + AgeDeath + onlyPTSD + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

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
save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_onlyPTSD_qsvs_withinDx.Rdata')


#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by dx
controlkeepIndex = which(rse_gene$onlyPTSD == "0")
rse_gene <- rse_gene[, controlkeepIndex]

#Get rid of dx columns
colIndex <-  !grepl("onlyPTSD", colnames(mod))
mod <- mod[controlkeepIndex, colIndex]
modQsva <- modQsva[controlkeepIndex,!grepl("onlyPTSD", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)


#Sex effects
withinDx_controls = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(withinDx_controls) = paste0(colnames(withinDx_controls), "_controlsonlyPTSD")





############# PTSD ##############


#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_withinDx.Rdata')


#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by dx
onlyptsdkeepIndex = which(rse_gene$onlyPTSD == "1")
rse_gene <- rse_gene[, onlyptsdkeepIndex]

#Get rid of dx columns
colIndex <-  !grepl("onlyPTSD", colnames(mod))
mod <- mod[onlyptsdkeepIndex, colIndex]
modQsva <- modQsva[onlyptsdkeepIndex,!grepl("onlyPTSD", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)


#Sex effects
withinDx_ptsd = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(withinDx_ptsd) = paste0(colnames(withinDx_ptsd), "_PTSDonlyPTSD")



###Stuff we want

###### GENE ########


#All columns
geneStatsall_withinDx = cbind(withinDx_controls, withinDx_ptsd)
geneStatsall_withinDx = cbind(geneStatsall_withinDx, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall_withinDx, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_withinDx_onlyPTSD.rda")

write.csv(geneStatsall_withinDx, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_withinDx_onlyPTSD.csv.gz"))



