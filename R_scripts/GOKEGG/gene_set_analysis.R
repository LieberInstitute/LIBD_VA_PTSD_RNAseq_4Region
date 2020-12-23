library(clusterProfiler)
library(org.Hs.eg.db)

## load data
load("../rdas/geneStats_DE_qSVA_mergedRegions_threeGroup.rda",verbose=TRUE)

## overall
p005 = dat[,grep("P.Value", colnames(dat))] < 0.005
colnames(p005) = paste0(ss(colnames(p005), "\\."), "_", ss(colnames(p005), "_",2))

## by dir
dirIndex = grepl("P.Value", colnames(dat)) & !grepl("ANOVA", colnames(dat))
p005.up = dat[,dirIndex] < 0.005 & dat[,grep("logFC", colnames(dat))] > 0
colnames(p005.up)= paste0(colnames(p005.up), "_UP")
p005.down = dat[,dirIndex] < 0.005 & dat[,grep("logFC", colnames(dat))] < 0
colnames(p005.down)= paste0(colnames(p005.down), "_DOWN")

## checks
colSums(p005.up) + colSums(p005.down)
colSums(p005)

## make matrix
geneInclMat= cbind(p005, p005.up, p005.down)
geneList = apply(geneInclMat, 2, function(x) {
	o = dat$EntrezID[x]
	as.character(o[!is.na(o)])
})

# and universe of expressed genes
universe = as.character(dat$EntrezID[!is.na(dat$EntrezID)])

## run GO and KEGG
go <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
kegg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)
## save
save(go, kegg, file = "../rdas/geneSet_threeGroups_qSVA.rda")

## check adjust p-values
goCheck <- compareCluster(geneList[1:2], fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
				
#############################
####### check results #######
load("../rdas/geneSet_threeGroups_qSVA.rda")

