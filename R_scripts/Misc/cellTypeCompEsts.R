#cell type composition estimation
library(BRETIGEA)
library(SummarizedExperiment)
library(jaffelab)
library(recount)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

#####################################
#Amyg
#####################################

###########
# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

controls <- rownames(colData(rse_gene))[rse_gene$Group == "Control"]
mdd <- rownames(colData(rse_gene))[rse_gene$Group == "MDD"]
ptsd <- rownames(colData(rse_gene))[rse_gene$Group == "PTSD"]


geneExprs = log2(getRPKM(rse_gene,"Length")+1)
rownames(geneExprs) <-rowData(rse_gene)$Symbol

geneExprs_cont <- geneExprs[,colnames(geneExprs) %in% controls]
geneExprs_MDD <- geneExprs[,colnames(geneExprs) %in% mdd]
geneExprs_PTSD <- geneExprs[,colnames(geneExprs) %in% ptsd]

log_prop_cont <- brainCells(geneExprs_cont, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)
log_prop_MDD <- brainCells(geneExprs_MDD, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)
log_prop_PTSD <- brainCells(geneExprs_PTSD, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)

apply(log_prop_PTSD, 2, function(x) t.test(x,log_prop_cont))
apply(log_prop_PTSD, 2, function(x) t.test(x,log_prop_MDD))
apply(log_prop_MDD, 2, function(x) t.test(x,log_prop_cont))


##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,3,4,5,6,20,2,7:19,21:40)]
head(modQsva)

## clean expression
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

geneExprsClean_cont <- geneExprsClean[,colnames(geneExprsClean) %in% controls]
geneExprsClean_MDD <- geneExprsClean[,colnames(geneExprsClean) %in% mdd]
geneExprsClean_PTSD <- geneExprsClean[,colnames(geneExprsClean) %in% ptsd]

clog_prop_cont <- brainCells(geneExprsClean_cont, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)
clog_prop_MDD <- brainCells(geneExprsClean_MDD, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)
clog_prop_PTSD <- brainCells(geneExprsClean_PTSD, nMarker = 50, species = "combined",celltypes = c("ast", "end", "mic", "neu", "oli", "opc"),method = "SVD", scale = TRUE)

apply(clog_prop_PTSD, 2, function(x) t.test(x,clog_prop_cont))
apply(clog_prop_PTSD, 2, function(x) t.test(x,clog_prop_MDD))
apply(clog_prop_MDD, 2, function(x) t.test(x,clog_prop_cont))

#Doesnt seem to be what I want
#let's try CIBERSORTx using Li 2018 expression data
#cell type composition estimation
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(readr)
library(devtools)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

#####################################
#Amyg
#####################################

###########
# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

controls <- rownames(colData(rse_gene))[rse_gene$Group == "Control"]
ptsd <- rownames(colData(rse_gene))[rse_gene$Group == "PTSD"]

load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,3,4,5,6,20,2,7:19,21:40)]
head(modQsva)

## clean expression
geneExprs = getRPKM(rse_gene,"Length")
geneExprs = cleaningY(geneExprs, modQsva, P=6)

rownames(geneExprs) <- rowData(rse_gene)$Symbol

geneExprs_cont <- geneExprs[,colnames(geneExprs) %in% controls]
geneExprs_PTSD <- geneExprs[,colnames(geneExprs) %in% ptsd]

write.csv(geneExprs_cont,file="csvs/geneExprs_cleaned_Amyg_Controls.csv")
write.csv(geneExprs_PTSD,file="csvs/geneExprs_cleaned_Amyg_PTSD.csv")

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/Data/Li_2018/Sestan.adultHumanNuclei.Psychencode.Rdata",verbose=TRUE)
#the meta2 file in the above indicates which nuclei are left after QC filtering. ill just match those nuclei numbers to the raw counts matrix to subset. 

umi.raw.subset <- umi.raw[,which(colnames(umi.raw) %in% rownames(meta2))]
#get rownames correctly formatted as well
rownames(umi.raw.subset) <- sub(".*\\|","",rownames(umi.raw.subset))

#drop any genes with no expression across all cells/nuclei
min(umi.raw.subset)
sum(is.na(umi.raw.subset))
totals <- rowSums(umi.raw.subset)
min(totals)
#appears no genes have zero counts across all nuclei (any genes that did were probably dropped in Li 2018's QC), but do below to check
umi.raw.subset_notzero <- umi.raw.subset[rowSums(umi.raw.subset) > 0,]
dim(umi.raw.subset_notzero)
dim(umi.raw.subset)

write.csv(umi.raw.subset_notzero, "csvs/Li_2018_countsmatrix.csv")








