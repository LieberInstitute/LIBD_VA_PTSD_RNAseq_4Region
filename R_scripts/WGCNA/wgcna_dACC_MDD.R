#######################################################################
#dACC
#MDD

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

## multithread
allowWGCNAThreads(8)

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

###########
# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

##########
## model #
##########

load("rdas/PTSD_qsvs.Rdata",verbose=TRUE)

##ensure intercept, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
modQsva = modQsva[,c(1,2,4,5,3,6:39)]

#filter for dACC
#gene
keepIndex = which(rse_gene$Region == "dACC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=4)


#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,networkType = "signed", verbose = 5)

cat(sftthresh1$powerEstimate)

save(sftthresh1, file = "rdas/WGCNA/power_object_dACC_MDD.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "rdas/WGCNA/wgcna_signed_TOM_dACC_MDD")
fNames = rownames(geneExprs)
save(net, fNames, file = "rdas/WGCNA/constructed_network_signed_bicor_dACC_MDD.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()
