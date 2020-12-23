#####Stratify by drugs of abuse analysis#######

###opioids_tox


library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

#load demographics data
demo <- read.csv('csvs/rse_gene_PTSD_cohort_Feb2019.csv',header=TRUE)
demo$BrNum <- demo$brnum
rownames(demo) <- demo$BrNum


#merge any_drugs column of demo with rse_object
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,84, drop=FALSE])


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


## also add snpPCs based on association with group
mod = model.matrix(~opioids_tox + Group + Sex + AgeDeath + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

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

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_qsvs_opioidstox.Rdata')



#filter for BasoAmyg 
#gene
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

#drugs vs no drugs
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene) = paste0(colnames(sigGene), "_opioidstox")






###Stuff we want

###### GENE ########


#All columns
geneStatsall = cbind(sigGene, rowData(rse_gene))
## write out qSVA-based stats
geneStats_BasoAmygall_opioidstox = geneStatsall
save(geneStats_BasoAmygall_opioidstox, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_opioidstox.rda")

write.csv(geneStats_BasoAmygall_opioidstox, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_opioidstox.csv.gz"))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()