###
library(jaffelab)
library(GenomicRanges)
library(limma)
library(SummarizedExperiment)

## write out annotation
load("../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")
geneMap = as.data.frame(rowRanges(rse_gene))
geneMap$gencodeTx = geneMap$Class = NULL
colnames(geneMap)[1] = "chr_hg38"
colnames(geneMap)[2:3] = paste0(colnames(geneMap)[2:3], "_hg38")
write.csv(geneMap, file = "../csvs/geneMap_gencodeV25.csv")

## pick files
fileList = list.files("../csvs", pattern = "DE_qSVA", recur=TRUE, full = TRUE)
names(fileList) = ss(fileList, "/", 3)
## filter
fileList = fileList[c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg")]

datList = lapply(fileList, read.csv, as.is=TRUE, row.names=1)
sapply(datList, dim)

## drop annotation info
dat = do.call("cbind", lapply(datList, function(x) x[,1:19]))

dat$Symbol = datList$DLPFC$Symbol
dat = dat[,c(ncol(dat), 1:(ncol(dat)-1))]
dat = cbind(dat, datList$DLPFC[,c("ensemblID", "gene_type", "EntrezID", "NumTx")]) # more annotation
save(dat, file = "../rdas/geneStats_DE_qSVA_mergedRegions_threeGroup.rda")
write.csv(dat, file = gzfile("../csvs/geneStats_DE_qSVA_mergedRegions_threeGroup.csv.gz"))

## any signif
anySig = dat[rowSums(dat[,grep("P.Value", colnames(dat))] < 0.001) > 0,]
dim(anySig)
save(anySig, file = "../rdas/geneStats_DE_qSVA_mergedRegions_threeGroup_p001any.rda")
write.csv(anySig, file = gzfile("../csvs/geneStats_DE_qSVA_mergedRegions_threeGroup_p001any.csv.gz"))

## how many sig?
colSums(dat[,grep("P.Value", colnames(dat))] < 0.001)
colSums(dat[,grep("P.Value", colnames(dat))] < 0.005)
colSums(dat[,grep("P.Value", colnames(dat))] < 0.01)
colSums(dat[,grep("adj.P.Val", colnames(dat))] < 0.05)
colSums(dat[,grep("adj.P.Val", colnames(dat))] < 0.1)
colSums(dat[,grep("adj.P.Val", colnames(dat))] < 0.2)

## venns
#Code adjusted by BKB on 4 March 2019 to account for fifth analysis (PTSDvsMDD) and to add color to the venn diagrams
#Issue: grep not discerning between "_PTSD" and "_PTSDvsMDD" so changed code for grep as needed

pdf("../pdf/venn_diagrams_genes_by_dx_and_region.pdf")
par(mar=c(1,1,1,1), cex.axis=1.5,cex.lab=1.5,cex.main = 1.5)
## within region, across dx
vennDiagram(vennCounts((dat[,grep("adj.P.Val", colnames(dat))] < 0.1)[,1:5]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Baso Amyg: q < 0.1", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val", colnames(dat))] < 0.1)[,6:10]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "dACC: q < 0.1", circle.col = c("red","blue","green","black","orange"))
	vennDiagram(vennCounts((dat[,grep("adj.P.Val", colnames(dat))] < 0.1)[,11:15]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "DLPFC: q < 0.1", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val", colnames(dat))] < 0.1)[,16:20]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Medial Amyg: q < 0.1", circle.col = c("red","blue","green","black","orange"))
		
## not much overlap## within region, across dx
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.001)[,1:5]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Baso Amyg: p < 0.001", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.001)[,6:10]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "dACC: p < 0.001", circle.col = c("red","blue","green","black","orange"))
	vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.001)[,11:15]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "DLPFC: p < 0.001", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.001)[,16:20]),
	names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Medial Amyg: p < 0.001", circle.col = c("red","blue","green","black","orange"))
## not much overlap

## more liberal; within region, across dx
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.01)[,1:5]),
		names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Baso Amyg: p < 0.01", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.01)[,6:10]),
		names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "dACC: p < 0.01", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.01)[,11:15]),
		names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "DLPFC: p < 0.01", circle.col = c("red","blue","green","black","orange"))
vennDiagram(vennCounts((dat[,grep("P.Value", colnames(dat))] < 0.01)[,16:20]),
		names = c("PTSD", "MDD", "ANOVA", "onlyPTSD", "PTSDvsMDD"),main = "Medial Amyg: p < 0.01", circle.col = c("red","blue","green","black","orange"))
## more mdd and ptsd overlap

## across region, w/in dx
vennDiagram(vennCounts((dat[,grep("adj.P.Val_MDD", colnames(dat))] < 0.1)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "MDD: q < 0.1", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val_PTSD$", colnames(dat))] < 0.1)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSD: q < 0.1", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val_onlyPTSD", colnames(dat))] < 0.1)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "onlyPTSD: q < 0.1", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val_ANOVA", colnames(dat))] < 0.1)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "ANOVA: q < 0.1", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("adj.P.Val_PTSDvsMDD", colnames(dat))] < 0.1)),
        names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSDvsMDD: q < 0.1", circle.col = c("red","blue","green","black"))


## across region, w/in dx
vennDiagram(vennCounts((dat[,grep("P.Value_MDD", colnames(dat))] < 0.001)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "MDD: p < 0.001", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_PTSD$", colnames(dat))] < 0.001)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSD: p < 0.001", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_onlyPTSD", colnames(dat))] < 0.001)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "ANOVA: p < 0.001", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_ANOVA", colnames(dat))] < 0.001)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "OnlyPTSD: p < 0.001", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_PTSDvsMDD", colnames(dat))] < 0.001)),
        names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSDvsMDD: p < 0.001", circle.col = c("red","blue","green","black"))
## not much overlap


## across region, w/in dx
vennDiagram(vennCounts((dat[,grep("P.Value_MDD", colnames(dat))] < 0.01)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "MDD: p < 0.01", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_PTSD$", colnames(dat))] < 0.01)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSD: p < 0.01", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_onlyPTSD", colnames(dat))] < 0.01)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "onlyPTSD: p < 0.01", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_ANOVA", colnames(dat))] < 0.01)),
	names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "ANOVA: p < 0.01", circle.col = c("red","blue","green","black"))
vennDiagram(vennCounts((dat[,grep("P.Value_PTSDvsMDD", colnames(dat))] < 0.01)),
        names = c("BasoAmyg", "dACC", "DLPFC", "MedialAmyg"),main = "PTSDvsMDD: p < 0.01", circle.col = c("red","blue","green","black"))
dev.off()
## not much overlap

#####################
## gene ontology ####

## overall
p005 = dat[,grep("P.Value", colnames(dat))] < 0.005
colnames(p005) = paste0(ss(colnames(p005), "\\."), "_", ss(colnames(p005), "_",2))

## by dir
dirIndex = grepl("P.Value", colnames(dat)) & !grepl("ANOVA", colnames(dat))
p005.up = dat[,dirIndex] < 0.005 & dat[,grep("logFC", colnames(dat))] > 0
colnames(p005.up)= paste0(colnames(p005.up), "_UP")
p005.down = dat[,dirIndex] < 0.005 & dat[,grep("logFC", colnames(dat))] < 0
colnames(p005.down)= paste0(colnames(p005.down), "_DOWN")

colSums(p005.up) + colSums(p005.down)
colSums(p005)

geneInclMat= cbind(p005, p005.up, p005.down)
geneList = apply(geneInclMat, 2, function(x) {
	o = dat$EntrezID[x]
	as.character(o[!is.na(o)])
})

universe = as.character(dat$EntrezID[!is.na(dat$EntrezID)])

go <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
kegg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe,  pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)
save(go, kegg, file = "../rdas/geneSet_threeGroups_qSVA.rda")
