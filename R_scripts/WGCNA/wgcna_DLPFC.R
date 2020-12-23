#######################################################################
#DLPFC
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
load('Data/rdas/General/RSEs/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

###########
# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

##########
## model #
##########

load("Data/rdas/General/ModelMatrices/Main/EachRegion/PTSD_qsvs.Rdata",verbose=TRUE)

##ensure intercept, group, agedeath, and sex are up first for cleaning
colnames(modQsva)

#filter for DLPFC 
#gene
keepIndex = which(rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)


#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,networkType = "signed", verbose = 5)

cat(sftthresh1$powerEstimate)

save(sftthresh1, file = "Data/rdas/WGCNA/DLPFC/power_object_DLPFC.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "Data/rdas/WGCNA/DLPFC/wgcna_signed_TOM_DLPFC")
fNames = rownames(geneExprs)
save(net, fNames, file = "Data/rdas/WGCNA/DLPFC/constructed_network_signed_bicor_DLPFC.rda")


###########################
### test ME associations ##
###########################
coefAdj_MDD = t(apply(net$MEs, 2, function(y)
  summary(lm(y~ modQsva[ ,1:5] -1))$coef[2,]))
coefAdj_PTSD = t(apply(net$MEs, 2, function(y)
  summary(lm(y~ modQsva[ ,1:5] -1))$coef[3,]))
coefAdj_onlyPTSD = t(apply(net$MEs, 2, function(y)
  summary(lm(y~ modQsva[ ,c(1,3:5)] -1))$coef[2,]))

save(coefAdj_MDD,coefAdj_PTSD,coefAdj_onlyPTSD, file="Data/rdas/WGCNA/DLPFC/MEvsDx_DLPFC.rda")

#GO using clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
				  
#make matrix
# split genes into modules
moduleGeneList_adj = split(rowData(rse_gene)$EntrezID, net$colors)
moduleGeneList_adj = lapply(moduleGeneList_adj, function(x) x[!is.na(x)])
names(moduleGeneList_adj) = paste0("DLPFC_", names(moduleGeneList_adj))

#universe of expressed genes
universe <- as.character(rowData(rse_gene)$EntrezID[!is.na(rowData(rse_gene)$EntrezID)])

## run GO 
go_modules <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)

save(go_modules, file="Data/rdas/WGCNA/DLPFC/wgcna_DLPFC_GO_clusterProfiler.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()

q()
