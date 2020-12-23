library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)
library(lmerTest)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_PTSD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

###########
# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

##########
## model #
##########
load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,3,4,5,6,20,2,7:19,21:40)]
head(modQsva)

### Dx association ###
brnums<-rse_gene$BrNum

coefAdj = t(apply(net$MEs, 2, function(y)
	summary(lmerTest::lmer(y~ modQsva[ ,1:6] - 1 + (1|brnums)))$coef[2,]))

save(coefAdj, file="rdas/WGCNA/MEvsDx_Amyg_PTSD.rda")

coefAdj[coefAdj[,5] < 0.1,]

#GO using clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
				  
#make matrix
# split genes into modules
moduleGeneList_adj = split(rowData(rse_gene)$EntrezID, net$colors)
moduleGeneList_adj = lapply(moduleGeneList_adj, function(x) x[!is.na(x)])

#universe of expressed genes
universe <- as.character(rowData(rse_gene)$EntrezID[!is.na(rowData(rse_gene)$EntrezID)])

## run GO 
go_modules_amyg <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1, readable= TRUE)

save(go_modules_amyg, file="rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_PTSD.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()
