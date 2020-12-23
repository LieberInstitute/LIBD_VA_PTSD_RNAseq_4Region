#DLPFC
#antidep based on tox
#PTSD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_DE_qSVA_sensitivity_antidepressants_DLPFC_threeGroup.rda")
geneStats_DLPFC <- as.data.frame(geneStats_DLPFC)
DLPFC_PTSD <- geneStats_DLPFC[,grep("_PTSD\\b", colnames(geneStats_DLPFC))]
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


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD_sensitivity_antidepressants.rda")

#DLPFC
#lifetime antidep
#PTSD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_DE_qSVA_sensitivity_lifetimeantidepressants_DLPFC_threeGroup.rda")
geneStats_DLPFC <- as.data.frame(geneStats_DLPFC)
DLPFC_PTSD <- geneStats_DLPFC[,grep("_PTSD\\b", colnames(geneStats_DLPFC))]
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


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD_sensitivity_lifetimeantidepressants.rda")

#DLPFC
#antidep based on tox
#MDD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_DE_qSVA_sensitivity_antidepressants_DLPFC_threeGroup.rda")
geneStats_DLPFC <- as.data.frame(geneStats_DLPFC)
DLPFC_MDD <- geneStats_DLPFC[,grep("_MDD", colnames(geneStats_DLPFC))]
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


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD_sensitivity_antidepressants.rda")

#DLPFC
#lifetime antidep
#MDD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/DLPFC/geneStats_DE_qSVA_sensitivity_lifetimeantidepressants_DLPFC_threeGroup.rda")
geneStats_DLPFC <- as.data.frame(geneStats_DLPFC)
DLPFC_MDD <- geneStats_DLPFC[,grep("_MDD", colnames(geneStats_DLPFC))]
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


save(go_DLPFC, kegg_DLPFC, file="rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD_sensitivity_lifetimeantidepressants.rda")

#dACC
#antidep based on tox
#MDD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_antidepressants_dACC_threeGroup.rda")
geneStats_dACC <- as.data.frame(geneStats_dACC)
dACC_MDD <- geneStats_dACC[,grep("_MDD", colnames(geneStats_dACC))]
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


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD_sensitivity_antidepressants.rda")

#dACC
#lifetime antidepressants
#MDD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_lifetimeantidepressants_dACC_threeGroup.rda")
geneStats_dACC <- as.data.frame(geneStats_dACC)
dACC_MDD <- geneStats_dACC[,grep("_MDD", colnames(geneStats_dACC))]
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


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD_sensitivity_lifetimeantidepressants.rda")

#dACC
#antidep based on tox
#antideptox effect
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_antidepressants_dACC_threeGroup_antideptox.rda")
geneStats_dACC <- as.data.frame(geneStats_dACC_antideptox)
dACC_antidep <- geneStats_dACC[,grep("_AntidepTox", colnames(geneStats_dACC))]
## overall
p005all = dACC_antidep[,grep("P.Value", colnames(dACC_antidep))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_antidep)) 
p005up = dACC_antidep[,dirIndex] < 0.005 & dACC_antidep[,grep("logFC", colnames(dACC_antidep))] > 0
p005down = dACC_antidep[,dirIndex] < 0.005 & dACC_antidep[,grep("logFC", colnames(dACC_antidep))] < 0

## checks
sum(p005up + p005down)
sum(p005all)

## make matrix
geneInclMat= cbind(p005all, p005up, p005down)
geneList = apply(geneInclMat, 2, function(x) {
        o = dACC_antidep$EntrezID[x]
        as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dACC_antidep$EntrezID[!is.na(dACC_antidep$EntrezID)])

## run GO and KEGG
go_dACC <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)
kegg_dACC <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_antideptox_sensitivity_antidepressants.rda")

#dACC
#eliminating all positive (or NA) for antidepressants_ssris (tox)
#MDD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_notoxantidepressants_dACC_threeGroup.rda")
geneStats_dACC <- as.data.frame(geneStats_dACC)
dACC_MDD <- geneStats_dACC[,grep("_MDD", colnames(geneStats_dACC))]

## overall
p005all = dACC_MDD[,grep("P.Value", colnames(dACC_MDD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_MDD)) 
p005up = dACC_MDD[,dirIndex] < 0.005 & dACC_MDD[,grep("logFC", colnames(dACC_MDD))] > 0
p005down = dACC_MDD[,dirIndex] < 0.005 & dACC_MDD[,grep("logFC", colnames(dACC_MDD))] < 0


## checks
sum(p005up + p005down)
#437
sum(p005all)
#437

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


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD_sensitivity_notoxantidepressants.rda")

#dACC
#eliminating all positive (or NA) for antidepressants_ssris (tox)
#PTSD
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_notoxantidepressants_dACC_threeGroup.rda")
geneStats_dACC <- as.data.frame(geneStats_dACC)
dACC_PTSD <- geneStats_dACC[,grep("_PTSD\\b", colnames(geneStats_dACC))]
## overall
p005all = dACC_PTSD[,grep("P.Value", colnames(dACC_PTSD))] < 0.005

## by dir
dirIndex = grepl("P.Value", colnames(dACC_PTSD)) 
p005up = dACC_PTSD[,dirIndex] < 0.005 & dACC_PTSD[,grep("logFC", colnames(dACC_PTSD))] > 0
p005down = dACC_PTSD[,dirIndex] < 0.005 & dACC_PTSD[,grep("logFC", colnames(dACC_PTSD))] < 0

## checks
sum(p005up + p005down)
#516
sum(p005all)
#516

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


save(go_dACC, kegg_dACC, file="rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD_sensitivity_notoxantidepressants.rda")


#compare to original model (not incorporating antidepressants (tox)
#DLPFC, PTSD
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#original 
load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")
geneStats_DLPFCall <- as.data.frame(geneStats_DLPFCall)
DLPFC_PTSD <- geneStats_DLPFCall[,grep("_PTSD\\b", colnames(geneStats_DLPFCall))]

#antidep (tox)
load("rdas/DLPFC/geneStats_DE_qSVA_sensitivity_antidepressants_DLPFC_threeGroup.rda")
geneStats_DLPFC_antidep <- as.data.frame(geneStats_DLPFC)
DLPFC_PTSD_antidep <- geneStats_DLPFC_antidep[,grep("_PTSD\\b", colnames(geneStats_DLPFC_antidep))]

#change blank to NA and get rid of NA
DLPFC_PTSD[DLPFC_PTSD == ""] <- NA
DLPFC_PTSD_antidep[DLPFC_PTSD_antidep == ""] <- NA

DLPFC_PTSD <- DLPFC_PTSD[!is.na(DLPFC_PTSD$Symbol_PTSD),]
DLPFC_PTSD_antidep <- DLPFC_PTSD_antidep[!is.na(DLPFC_PTSD_antidep$Symbol_PTSD),]

DLPFC_PTSD_sig <- DLPFC_PTSD[DLPFC_PTSD$P.Value_PTSD < 0.005,]
#358
DLPFC_PTSD_antidep_sig <- DLPFC_PTSD_antidep[DLPFC_PTSD_antidep$P.Value_PTSD < 0.005,]
#308

DLPFC_PTSD_sig_up <- DLPFC_PTSD_sig[DLPFC_PTSD_sig$logFC_PTSD > 0,]
#201
DLPFC_PTSD_sig_down <- DLPFC_PTSD_sig[DLPFC_PTSD_sig$logFC_PTSD < 0,]
#157

DLPFC_PTSD_antidep_sig_up <- DLPFC_PTSD_antidep_sig[DLPFC_PTSD_antidep_sig$logFC_PTSD > 0,]
#191
DLPFC_PTSD_antidep_sig_down <- DLPFC_PTSD_antidep_sig[DLPFC_PTSD_antidep_sig$logFC_PTSD < 0,]
#118

length(unique(DLPFC_PTSD_sig$Symbol_PTSD)[which(DLPFC_PTSD_sig$Symbol_PTSD %in% DLPFC_PTSD_antidep_sig$Symbol_PTSD)])
#194
length(unique(DLPFC_PTSD_sig_up$Symbol_PTSD)[which(DLPFC_PTSD_sig_up$Symbol_PTSD %in% DLPFC_PTSD_antidep_sig_up$Symbol_PTSD)])
#116
length(unique(DLPFC_PTSD_sig_down$Symbol_PTSD)[which(DLPFC_PTSD_sig_down$Symbol_PTSD %in% DLPFC_PTSD_antidep_sig_down$Symbol_PTSD)])
#78

#not overlapping
up_non <- unique(setdiff(DLPFC_PTSD_sig_up$Symbol_PTSD,DLPFC_PTSD_antidep_sig_up$Symbol_PTSD))
#84
down_non <- unique(setdiff(DLPFC_PTSD_sig_down$Symbol_PTSD,DLPFC_PTSD_antidep_sig_down$Symbol_PTSD))
#79

DLPFC_PTSD_sig2 <- DLPFC_PTSD_sig[which(DLPFC_PTSD_sig$Symbol_PTSD %in% DLPFC_PTSD_antidep_sig$Symbol_PTSD),]
DLPFC_PTSD_antidep_sig2 <- DLPFC_PTSD_antidep_sig[match(DLPFC_PTSD_sig2$Symbol_PTSD, DLPFC_PTSD_antidep_sig$Symbol_PTSD),]
colnames(DLPFC_PTSD_antidep_sig2) <-paste0(colnames(DLPFC_PTSD_antidep_sig2),"_antidep")

both_sig <- cbind(DLPFC_PTSD_sig2,DLPFC_PTSD_antidep_sig2)

both_sig[which(both_sig$logFC_PTSD > 0 & both_sig$logFC_PTSD_antidep < 0),]
#0
both_sig[which(both_sig$logFC_PTSD < 0 & both_sig$logFC_PTSD_antidep > 0),]
#0

#compare to original model (not incorporating antidep tox) 
#dACC, MDD
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#original 
load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
geneStats_dACCall <- as.data.frame(geneStats_dACCall)
dACC_MDD <- geneStats_dACCall[,grep("_MDD", colnames(geneStats_dACCall))]

#antidep (tox)
load("rdas/dACC/geneStats_DE_qSVA_sensitivity_antidepressants_dACC_threeGroup.rda")
geneStats_dACC_antidep <- as.data.frame(geneStats_dACC)
dACC_MDD_antidep <- geneStats_dACC_antidep[,grep("_MDD", colnames(geneStats_dACC_antidep))]

#change blank to NA and get rid of NA
dACC_MDD[dACC_MDD == ""] <- NA
dACC_MDD_antidep[dACC_MDD_antidep == ""] <- NA

dACC_MDD <- dACC_MDD[!is.na(dACC_MDD$Symbol_MDD),]
dACC_MDD_antidep <- dACC_MDD_antidep[!is.na(dACC_MDD_antidep$Symbol_MDD),]

dACC_MDD_sig <- dACC_MDD[dACC_MDD$P.Value_MDD < 0.005,]
#668
dACC_MDD_antidep_sig <- dACC_MDD_antidep[dACC_MDD_antidep$P.Value_MDD < 0.005,]
#410

dACC_MDD_sig_up <- dACC_MDD_sig[dACC_MDD_sig$logFC_MDD > 0,]
#223
dACC_MDD_sig_down <- dACC_MDD_sig[dACC_MDD_sig$logFC_MDD < 0,]
#445

dACC_MDD_antidep_sig_up <- dACC_MDD_antidep_sig[dACC_MDD_antidep_sig$logFC_MDD > 0,]
#182
dACC_MDD_antidep_sig_down <- dACC_MDD_antidep_sig[dACC_MDD_antidep_sig$logFC_MDD < 0,]
#228

length(unique(dACC_MDD_sig$Symbol_MDD)[which(dACC_MDD_sig$Symbol_MDD %in% dACC_MDD_antidep_sig$Symbol_MDD)])
#322
length(unique(dACC_MDD_sig_up$Symbol_MDD)[which(dACC_MDD_sig_up$Symbol_MDD %in% dACC_MDD_antidep_sig_up$Symbol_MDD)])
#138
length(unique(dACC_MDD_sig_down$Symbol_MDD)[which(dACC_MDD_sig_down$Symbol_MDD %in% dACC_MDD_antidep_sig_down$Symbol_MDD)])
#183

#not overlapping
up_non <- unique(setdiff(dACC_MDD_sig_up$Symbol_MDD,dACC_MDD_antidep_sig_up$Symbol_MDD))
#85
down_non <- unique(setdiff(dACC_MDD_sig_down$Symbol_MDD,dACC_MDD_antidep_sig_down$Symbol_MDD))
#262

dACC_MDD_sig2 <- dACC_MDD_sig[which(dACC_MDD_sig$Symbol_MDD %in% dACC_MDD_antidep_sig$Symbol_MDD),]
dACC_MDD_antidep_sig2 <- dACC_MDD_antidep_sig[match(dACC_MDD_sig2$Symbol_MDD, dACC_MDD_antidep_sig$Symbol_MDD),]
colnames(dACC_MDD_antidep_sig2) <-paste0(colnames(dACC_MDD_antidep_sig2),"_antidep")

both_sig <- cbind(dACC_MDD_sig2,dACC_MDD_antidep_sig2)

both_sig[which(both_sig$logFC_MDD > 0 & both_sig$logFC_MDD_antidep < 0),]
#0
both_sig[which(both_sig$logFC_MDD < 0 & both_sig$logFC_MDD_antidep > 0),]
#1
#ENSG00000277118, misc_RNA
