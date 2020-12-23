library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

#Read in Joel's data
joel <- read.csv("/dcl02/lieber/ajaffe/kleinman_PTSD/Kleinman_R01_PTSD_RNA_Samples_Demo.csv")
load("/dcl02/lieber/ajaffe/kleinman_PTSD/preprocessed_data/rse_gene_Kleinman_R01_PTSD_n225.Rdata")

#set row names and bind sample data with rse_gene object
rownames(joel) = joel$RNum
colData(rse_gene) = cbind(colData(rse_gene), joel[rse_gene$SAMPLE_ID,])

#assign unique name
joel_rse_gene <- rse_gene

#Read in our data and assign unique name
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
og_rse_gene <- rse_gene

#Keep regions of interest
ogkeepIndex = which(og_rse_gene$Region == "BasoAmyg")
og_rse_gene <- og_rse_gene[,ogkeepIndex]
joelkeepIndex = which(joel_rse_gene$Brain.Region == "Hippocampus")
joel_rse_gene <- joel_rse_gene[,joelkeepIndex]

#Only keep samples where we have both the HPC and the BLA
test <- which(og_rse_gene$BrNum %in% joel_rse_gene$BrNum)
both_og_rse_gene <- og_rse_gene[,test]
test2 <- which(joel_rse_gene$BrNum %in% both_og_rse_gene$BrNum)
both_joel_rse_gene <- joel_rse_gene[,test2]

#Check that the brain numbers are the same
both_joel_rse_gene$BrNum <- as.character(both_joel_rse_gene$BrNum)
unique(both_joel_rse_gene$BrNum)
unique(both_og_rse_gene$BrNum)
identical(sort(both_joel_rse_gene$BrNum), sort(both_og_rse_gene$BrNum))

#Subset into the different diagnoses
summary(as.factor(both_og_rse_gene$Group))
summary(both_joel_rse_gene$PrimaryDx)

ptsdkeepIndex <- which(both_og_rse_gene$Group == "PTSD")
ptsd_both_og_rse_gene <- both_og_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_og_rse_gene$Group == "MDD")
mdd_both_og_rse_gene <- both_og_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_og_rse_gene$Group == "Control")
cont_both_og_rse_gene <- both_og_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_joel_rse_gene$PrimaryDx == "PTSD")
ptsd_both_joel_rse_gene <- both_joel_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_joel_rse_gene$PrimaryDx == "MDD")
mdd_both_joel_rse_gene <- both_joel_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_joel_rse_gene$PrimaryDx == "Control")
cont_both_joel_rse_gene <- both_joel_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_joel_rse_gene)
dim(cont_both_og_rse_gene)
dim(ptsd_both_og_rse_gene)
dim(ptsd_both_joel_rse_gene)
dim(mdd_both_joel_rse_gene)
dim(mdd_both_og_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_og_rse_gene <- rowMeans(assay(ptsd_both_og_rse_gene))
means_mdd_both_og_rse_gene <- rowMeans(assay(mdd_both_og_rse_gene))
means_cont_both_og_rse_gene <- rowMeans(assay(cont_both_og_rse_gene))
means_cont_both_joel_rse_gene <- rowMeans(assay(cont_both_joel_rse_gene))
means_mdd_both_joel_rse_gene <- rowMeans(assay(mdd_both_joel_rse_gene))
means_ptsd_both_joel_rse_gene <- rowMeans(assay(ptsd_both_joel_rse_gene))

#Calculate Pearson correlations
cor(means_ptsd_both_joel_rse_gene,means_ptsd_both_og_rse_gene)
cor(means_mdd_both_joel_rse_gene,means_mdd_both_og_rse_gene)
cor(means_cont_both_joel_rse_gene,means_cont_both_og_rse_gene)

#subset to genes of interest, then check correlation? 
