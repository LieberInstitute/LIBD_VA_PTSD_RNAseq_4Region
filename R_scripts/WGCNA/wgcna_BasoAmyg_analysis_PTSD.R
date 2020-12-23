library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load network
load("rdas/WGCNA/constructed_network_signed_bicor_BasoAmyg_PTSD.rda", verbose = TRUE)
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
modQsva = modQsva[,c(1,3,4,5,2,6:39)]

#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"

### Dx association ###

coefAdj = t(apply(net$MEs, 2, function(y)
  summary(lm(y~ modQsva[ ,1:4]))$coef[2,]))

save(coefAdj, file="rdas/WGCNA/MEvsDx_BasoAmyg_PTSD.rda")

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
go_modules_basoamyg <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                                readable= TRUE)

save(go_modules_basoamyg, file="rdas/BasoAmyg/wgcna_BasoAmyg_GO_clusterProfiler_PTSD.rda")

#sig <- coefAdj[coefAdj[,4] < 0.001,]
#go_df <- as.data.frame(go_modules_basoamyg)
#rownames(sig) <- gsub("ME","",rownames(sig))
#go_interest <- go_df[go_df$Cluster %in% rownames(sig),] #because those modules are ~significant

#go_interest_ordered <- go_interest[order(go_interest$p.adjust,decreasing=FALSE),]
#go_interest_ordered$Cluster<-as.character(paste0("ME",go_interest_ordered$Cluster))
#go_interest_ordered_split <- split(go_interest_ordered, go_interest_ordered$Cluster)
#basoamyg_desc<-sapply(go_interest_ordered_split, "[","Description")
#basoamyg_desc_top <- sapply(basoamyg_desc, "[",c(1:5))

#basoamyg_desc[["ME12.Description"]][1:10]

#both<-sapply(go_interest_ordered_split, "[",c("Description","p.adjust"))
#both_top <- sapply(both, "[",c(1:10))
#go_interest_ordered_split[["ME12"]][c("Description","p.adjust","geneID")][1:10,]



							