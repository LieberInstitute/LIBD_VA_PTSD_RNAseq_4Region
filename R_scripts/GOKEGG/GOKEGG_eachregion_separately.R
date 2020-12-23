###################################################
##Gene ontology and KEGG broken down by region using MDD

#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")
geneStats_BasoAmygall <- as.data.frame(geneStats_BasoAmygall)
BLA_MDD <- geneStats_BasoAmygall[,grep("_MDD", colnames(geneStats_BasoAmygall))]
## overall
p005all = BLA_MDD[,grep("P.Value", colnames(BLA_MDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(BLA_MDD)) 
p005up = BLA_MDD[,dirIndex] < 0.005 & BLA_MDD[,grep("logFC", colnames(BLA_MDD))] > 0
p005down = BLA_MDD[,dirIndex] < 0.005 & BLA_MDD[,grep("logFC", colnames(BLA_MDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = BLA_MDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(BLA_MDD$EntrezID[!is.na(BLA_MDD$EntrezID)])

## run GO and KEGG
go_basoamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_basoamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_basoamyg, kegg_basoamyg, file="rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_MDD.rda")

###################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda")
geneStats_MedialAmygall <- as.data.frame(geneStats_MedialAmygall)
MeA_MDD <- geneStats_MedialAmygall[,grep("_MDD", colnames(geneStats_MedialAmygall))]
## overall
p005all = MeA_MDD[,grep("P.Value", colnames(MeA_MDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(MeA_MDD)) 
p005up = MeA_MDD[,dirIndex] < 0.005 & MeA_MDD[,grep("logFC", colnames(MeA_MDD))] > 0
p005down = MeA_MDD[,dirIndex] < 0.005 & MeA_MDD[,grep("logFC", colnames(MeA_MDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = MeA_MDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(MeA_MDD$EntrezID[!is.na(MeA_MDD$EntrezID)])

## run GO and KEGG
go_medialamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_medialamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_medialamyg, kegg_medialamyg, file="rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_MDD.rda")

###################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
geneStats_dACCall <- as.data.frame(geneStats_dACCall)
dACC_MDD <- geneStats_dACCall[,grep("_MDD", colnames(geneStats_dACCall))]
## overall
p005all = dACC_MDD[,grep("P.Value", colnames(dACC_MDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_MDD)) 
p005up = dACC_MDD[,dirIndex] < 0.005 & dACC_MDD[,grep("logFC", colnames(dACC_MDD))] > 0
p005down = dACC_MDD[,dirIndex] < 0.005 & dACC_MDD[,grep("logFC", colnames(dACC_MDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = dACC_MDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dACC_MDD$EntrezID[!is.na(dACC_MDD$EntrezID)])

## run GO and KEGG
go_dACC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_dACC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD.rda")

###################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")
geneStats_DLPFCall <- as.data.frame(geneStats_DLPFCall)
DLPFC_MDD <- geneStats_DLPFCall[,grep("_MDD", colnames(geneStats_DLPFCall))]
## overall
p005all = DLPFC_MDD[,grep("P.Value", colnames(DLPFC_MDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(DLPFC_MDD)) 
p005up = DLPFC_MDD[,dirIndex] < 0.005 & DLPFC_MDD[,grep("logFC", colnames(DLPFC_MDD))] > 0
p005down = DLPFC_MDD[,dirIndex] < 0.005 & DLPFC_MDD[,grep("logFC", colnames(DLPFC_MDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = DLPFC_MDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(DLPFC_MDD$EntrezID[!is.na(DLPFC_MDD$EntrezID)])

## run GO and KEGG
go_DLPFC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_DLPFC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD.rda")

#########################################################################
#PTSD

#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")
geneStats_BasoAmygall <- as.data.frame(geneStats_BasoAmygall)
BLA_PTSD <- geneStats_BasoAmygall[,grep("_PTSD\\b", colnames(geneStats_BasoAmygall))]
## overall
p005all = BLA_PTSD[,grep("P.Value", colnames(BLA_PTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(BLA_PTSD)) 
p005up = BLA_PTSD[,dirIndex] < 0.005 & BLA_PTSD[,grep("logFC", colnames(BLA_PTSD))] > 0
p005down = BLA_PTSD[,dirIndex] < 0.005 & BLA_PTSD[,grep("logFC", colnames(BLA_PTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = BLA_PTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(BLA_PTSD$EntrezID[!is.na(BLA_PTSD$EntrezID)])

## run GO and KEGG
go_basoamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_basoamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_basoamyg, kegg_basoamyg, file="rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSD.rda")

###################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda")
geneStats_MedialAmygall <- as.data.frame(geneStats_MedialAmygall)
MeA_PTSD <- geneStats_MedialAmygall[,grep("_PTSD\\b", colnames(geneStats_MedialAmygall))]
## overall
p005all = MeA_PTSD[,grep("P.Value", colnames(MeA_PTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(MeA_PTSD)) 
p005up = MeA_PTSD[,dirIndex] < 0.005 & MeA_PTSD[,grep("logFC", colnames(MeA_PTSD))] > 0
p005down = MeA_PTSD[,dirIndex] < 0.005 & MeA_PTSD[,grep("logFC", colnames(MeA_PTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = MeA_PTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(MeA_PTSD$EntrezID[!is.na(MeA_PTSD$EntrezID)])

## run GO and KEGG
go_medialamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_medialamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_medialamyg, kegg_medialamyg, file="rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSD.rda")

###################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
geneStats_dACCall <- as.data.frame(geneStats_dACCall)
dACC_PTSD <- geneStats_dACCall[,grep("_PTSD\\b", colnames(geneStats_dACCall))]
## overall
p005all = dACC_PTSD[,grep("P.Value", colnames(dACC_PTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_PTSD)) 
p005up = dACC_PTSD[,dirIndex] < 0.005 & dACC_PTSD[,grep("logFC", colnames(dACC_PTSD))] > 0
p005down = dACC_PTSD[,dirIndex] < 0.005 & dACC_PTSD[,grep("logFC", colnames(dACC_PTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = dACC_PTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dACC_PTSD$EntrezID[!is.na(dACC_PTSD$EntrezID)])

## run GO and KEGG
go_dACC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_dACC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD.rda")

###################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")
geneStats_DLPFCall <- as.data.frame(geneStats_DLPFCall)
DLPFC_PTSD <- geneStats_DLPFCall[,grep("_PTSD\\b", colnames(geneStats_DLPFCall))]
## overall
p005all = DLPFC_PTSD[,grep("P.Value", colnames(DLPFC_PTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(DLPFC_PTSD)) 
p005up = DLPFC_PTSD[,dirIndex] < 0.005 & DLPFC_PTSD[,grep("logFC", colnames(DLPFC_PTSD))] > 0
p005down = DLPFC_PTSD[,dirIndex] < 0.005 & DLPFC_PTSD[,grep("logFC", colnames(DLPFC_PTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = DLPFC_PTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(DLPFC_PTSD$EntrezID[!is.na(DLPFC_PTSD$EntrezID)])

## run GO and KEGG
go_DLPFC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_DLPFC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD.rda")

####################################################################################
#onlyPTSD

#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")
geneStats_BasoAmygall <- as.data.frame(geneStats_BasoAmygall)
BLA_onlyPTSD <- geneStats_BasoAmygall[,grep("_onlyPTSD", colnames(geneStats_BasoAmygall))]
## overall
p005all = BLA_onlyPTSD[,grep("P.Value", colnames(BLA_onlyPTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(BLA_onlyPTSD)) 
p005up = BLA_onlyPTSD[,dirIndex] < 0.005 & BLA_onlyPTSD[,grep("logFC", colnames(BLA_onlyPTSD))] > 0
p005down = BLA_onlyPTSD[,dirIndex] < 0.005 & BLA_onlyPTSD[,grep("logFC", colnames(BLA_onlyPTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = BLA_onlyPTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(BLA_onlyPTSD$EntrezID[!is.na(BLA_onlyPTSD$EntrezID)])

## run GO and KEGG
go_basoamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_basoamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_basoamyg, kegg_basoamyg, file="rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_onlyPTSD.rda")

###################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda")
geneStats_MedialAmygall <- as.data.frame(geneStats_MedialAmygall)
MeA_onlyPTSD <- geneStats_MedialAmygall[,grep("_onlyPTSD", colnames(geneStats_MedialAmygall))]
## overall
p005all = MeA_onlyPTSD[,grep("P.Value", colnames(MeA_onlyPTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(MeA_onlyPTSD)) 
p005up = MeA_onlyPTSD[,dirIndex] < 0.005 & MeA_onlyPTSD[,grep("logFC", colnames(MeA_onlyPTSD))] > 0
p005down = MeA_onlyPTSD[,dirIndex] < 0.005 & MeA_onlyPTSD[,grep("logFC", colnames(MeA_onlyPTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = MeA_onlyPTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(MeA_onlyPTSD$EntrezID[!is.na(MeA_onlyPTSD$EntrezID)])

## run GO and KEGG
go_medialamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_medialamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_medialamyg, kegg_medialamyg, file="rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_onlyPTSD.rda")

###################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
geneStats_dACCall <- as.data.frame(geneStats_dACCall)
dACC_onlyPTSD <- geneStats_dACCall[,grep("_onlyPTSD", colnames(geneStats_dACCall))]
## overall
p005all = dACC_onlyPTSD[,grep("P.Value", colnames(dACC_onlyPTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_onlyPTSD)) 
p005up = dACC_onlyPTSD[,dirIndex] < 0.005 & dACC_onlyPTSD[,grep("logFC", colnames(dACC_onlyPTSD))] > 0
p005down = dACC_onlyPTSD[,dirIndex] < 0.005 & dACC_onlyPTSD[,grep("logFC", colnames(dACC_onlyPTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = dACC_onlyPTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dACC_onlyPTSD$EntrezID[!is.na(dACC_onlyPTSD$EntrezID)])

## run GO and KEGG
go_dACC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_dACC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_onlyPTSD.rda")

###################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")
geneStats_DLPFCall <- as.data.frame(geneStats_DLPFCall)
DLPFC_onlyPTSD <- geneStats_DLPFCall[,grep("_onlyPTSD", colnames(geneStats_DLPFCall))]
## overall
p005all = DLPFC_onlyPTSD[,grep("P.Value", colnames(DLPFC_onlyPTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(DLPFC_onlyPTSD)) 
p005up = DLPFC_onlyPTSD[,dirIndex] < 0.005 & DLPFC_onlyPTSD[,grep("logFC", colnames(DLPFC_onlyPTSD))] > 0
p005down = DLPFC_onlyPTSD[,dirIndex] < 0.005 & DLPFC_onlyPTSD[,grep("logFC", colnames(DLPFC_onlyPTSD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = DLPFC_onlyPTSD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(DLPFC_onlyPTSD$EntrezID[!is.na(DLPFC_onlyPTSD$EntrezID)])

## run GO and KEGG
go_DLPFC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_DLPFC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_onlyPTSD.rda")


###################################################
##Gene ontology and KEGG broken down by region using PTSDvsMDD

#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")
geneStats_BasoAmygall <- as.data.frame(geneStats_BasoAmygall)
BLA_PTSDvsMDD <- geneStats_BasoAmygall[,grep("_PTSDvsMDD", colnames(geneStats_BasoAmygall))]
## overall
p005all = BLA_PTSDvsMDD[,grep("P.Value", colnames(BLA_PTSDvsMDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(BLA_PTSDvsMDD)) 
p005up = BLA_PTSDvsMDD[,dirIndex] < 0.005 & BLA_PTSDvsMDD[,grep("logFC", colnames(BLA_PTSDvsMDD))] > 0
p005down = BLA_PTSDvsMDD[,dirIndex] < 0.005 & BLA_PTSDvsMDD[,grep("logFC", colnames(BLA_PTSDvsMDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = BLA_PTSDvsMDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(BLA_PTSDvsMDD$EntrezID[!is.na(BLA_PTSDvsMDD$EntrezID)])

## run GO and KEGG
go_basoamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_basoamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_basoamyg, kegg_basoamyg, file="rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSDvsMDD.rda")

###################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda")
geneStats_MedialAmygall <- as.data.frame(geneStats_MedialAmygall)
MeA_PTSDvsMDD <- geneStats_MedialAmygall[,grep("_PTSDvsMDD", colnames(geneStats_MedialAmygall))]
## overall
p005all = MeA_PTSDvsMDD[,grep("P.Value", colnames(MeA_PTSDvsMDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(MeA_PTSDvsMDD)) 
p005up = MeA_PTSDvsMDD[,dirIndex] < 0.005 & MeA_PTSDvsMDD[,grep("logFC", colnames(MeA_PTSDvsMDD))] > 0
p005down = MeA_PTSDvsMDD[,dirIndex] < 0.005 & MeA_PTSDvsMDD[,grep("logFC", colnames(MeA_PTSDvsMDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = MeA_PTSDvsMDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(MeA_PTSDvsMDD$EntrezID[!is.na(MeA_PTSDvsMDD$EntrezID)])

## run GO and KEGG
go_medialamyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_medialamyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_medialamyg, kegg_medialamyg, file="rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSDvsMDD.rda")

###################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
geneStats_dACCall <- as.data.frame(geneStats_dACCall)
dACC_PTSDvsMDD <- geneStats_dACCall[,grep("_PTSDvsMDD", colnames(geneStats_dACCall))]
## overall
p005all = dACC_PTSDvsMDD[,grep("P.Value", colnames(dACC_PTSDvsMDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_PTSDvsMDD)) 
p005up = dACC_PTSDvsMDD[,dirIndex] < 0.005 & dACC_PTSDvsMDD[,grep("logFC", colnames(dACC_PTSDvsMDD))] > 0
p005down = dACC_PTSDvsMDD[,dirIndex] < 0.005 & dACC_PTSDvsMDD[,grep("logFC", colnames(dACC_PTSDvsMDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = dACC_PTSDvsMDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dACC_PTSDvsMDD$EntrezID[!is.na(dACC_PTSDvsMDD$EntrezID)])

## run GO and KEGG
go_dACC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_dACC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSDvsMDD.rda")

###################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")
geneStats_DLPFCall <- as.data.frame(geneStats_DLPFCall)
DLPFC_PTSDvsMDD <- geneStats_DLPFCall[,grep("_PTSDvsMDD", colnames(geneStats_DLPFCall))]
## overall
p005all = DLPFC_PTSDvsMDD[,grep("P.Value", colnames(DLPFC_PTSDvsMDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(DLPFC_PTSDvsMDD)) 
p005up = DLPFC_PTSDvsMDD[,dirIndex] < 0.005 & DLPFC_PTSDvsMDD[,grep("logFC", colnames(DLPFC_PTSDvsMDD))] > 0
p005down = DLPFC_PTSDvsMDD[,dirIndex] < 0.005 & DLPFC_PTSDvsMDD[,grep("logFC", colnames(DLPFC_PTSDvsMDD))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = DLPFC_PTSDvsMDD$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(DLPFC_PTSDvsMDD$EntrezID[!is.na(DLPFC_PTSDvsMDD$EntrezID)])

## run GO and KEGG
go_DLPFC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_DLPFC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSDvsMDD.rda")


