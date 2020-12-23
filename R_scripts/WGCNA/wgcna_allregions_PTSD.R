
#WGCNA
#across regions
#PTSD

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
load('rdas/PTSD_qsvs_Regionintxn.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
head(modQsva)

modQsva = modQsva[,c(1,3:8,22,24,26,2,9:21,23,25,27:46)]

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=10)


#########################
## get power
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), networkType = "signed", verbose = 5)
cat(sftthresh1$powerEstimate)
save(sftthresh1, file = "rdas/WGCNA/power_object_allregions_PTSD.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "rdas/WGCNA/wgcna_signed_TOM_allregions_PTSD")
fNames = rownames(geneExprs)
save(net, fNames, file = "rdas/WGCNA/constructed_network_signed_bicor_allregions_PTSD.rda")

