library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load network
load("rdas/WGCNA/constructed_network_signed_bicor_DLPFC_onlyPTSD.rda", verbose = TRUE)
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

load("rdas/PTSD_onlyPTSD_qsvs.Rdata",verbose=TRUE)

##ensure intercept, group, agedeath, and sex are up first for cleaning
colnames(modQsva)

#filter for DLPFC
#gene
keepIndex = which(rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

### Dx association ###

coefAdj = t(apply(net$MEs, 2, function(y)
  summary(lm(y~ modQsva[ ,1:4]))$coef[2,]))

save(coefAdj, file="rdas/WGCNA/MEvsDx_DLPFC_onlyPTSD.rda")

#Prelim GO using WGCNA package functionality
#library(GO.db)
#library(org.Hs.eg.db)

#universe <- as.character(rowData(rse_gene)$EntrezID[!is.na(rowData(rse_gene)$EntrezID)])
#GOenr = GOenrichmentAnalysis(net$colors, universe, organism = "human", nBestP = 10);
#tab = GOenr$bestPTerms[[4]]$enrichment
#sig <- coefAdj[coefAdj[,4] < 0.1,]
#rownames(sig) <- gsub("ME","",rownames(sig))
#tab_interest <- tab[tab$module == rownames(sig),] #because those modules are ~significant
#unique(tab_interest$termName)

#save(GOenr, file="rdas/DLPFC/wgcna_DLPFC_GO_onlyPTSD.rda")

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
go_modules_dlpfc <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)

save(go_modules_dlpfc, file="rdas/DLPFC/wgcna_DLPFC_GO_clusterProfiler_onlyPTSD.rda")


#sig <- coefAdj[coefAdj[,4] < 0.1,]
#go_df <- as.data.frame(go_modules_basoamyg)
#rownames(sig) <- gsub("ME","",rownames(sig))
#go_interest <- go_df[go_df$Cluster == rownames(sig),] #because those modules are ~significant
#go_interest_padjust <- go_interest[go_interest$p.adjust <0.1,]
