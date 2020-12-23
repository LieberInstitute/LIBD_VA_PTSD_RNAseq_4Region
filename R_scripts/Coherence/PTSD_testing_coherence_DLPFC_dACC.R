##################################
##################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)


#Read in our data and assign unique name
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


#Keep regions of interest
dACCkeepIndex = which(rse_gene$Region == "dACC")
dACC_rse_gene <- rse_gene[,dACCkeepIndex]
DLPFCkeepIndex = which(rse_gene$Region == "DLPFC")
DLPFC_rse_gene <- rse_gene[,DLPFCkeepIndex]

#Only keep samples where we have both the dACC and the DLPFC
test <- which(dACC_rse_gene$BrNum %in% DLPFC_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(DLPFC_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_DLPFC_rse_gene <- DLPFC_rse_gene[,test2]

#Check that the brain numbers are the same
both_DLPFC_rse_gene$BrNum <- as.character(both_DLPFC_rse_gene$BrNum)
unique(both_DLPFC_rse_gene$BrNum)
unique(both_dACC_rse_gene$BrNum)
identical(sort(both_DLPFC_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_DLPFC_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_DLPFC_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_DLPFC_rse_gene$Group == "PTSD")
ptsd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_DLPFC_rse_gene$Group == "MDD")
mdd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_DLPFC_rse_gene$Group == "Control")
cont_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_DLPFC_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_DLPFC_rse_gene)
dim(mdd_both_DLPFC_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(assay(ptsd_both_dACC_rse_gene))
means_mdd_both_dACC_rse_gene <- rowMeans(assay(mdd_both_dACC_rse_gene))
means_cont_both_dACC_rse_gene <- rowMeans(assay(cont_both_dACC_rse_gene))
means_cont_both_DLPFC_rse_gene <- rowMeans(assay(cont_both_DLPFC_rse_gene))
means_mdd_both_DLPFC_rse_gene <- rowMeans(assay(mdd_both_DLPFC_rse_gene))
means_ptsd_both_DLPFC_rse_gene <- rowMeans(assay(ptsd_both_DLPFC_rse_gene))

#Calculate Pearson correlations
cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dACC_rse_gene)

#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_DLPFC_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_DLPFC_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_DLPFC_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

#subset to genes of interest, then check correlation? 
#maybe look at exons, other features?

#ideas for genes to subset to
##BDNF, FKBP5, Mtor, Akt, Egr, NMDARs, AMPARs, CREB, NTRK2
###USF, Sp3/4
####mGluR2, 5HT2AR

#gencode IDs
plasticity <- c("ENSG00000176697.18","ENSG00000096060.14","ENSG00000198793.12","ENSG00000120738.7","ENSG00000164082.14","ENSG00000102468.10","ENSG00000176884.14","ENSG00000148053.15","ENSG00000118260.14")

means_ptsd_both_dACC_rse_gene <- means_ptsd_both_dACC_rse_gene[which(names(means_ptsd_both_dACC_rse_gene) %in% plasticity)] 
means_mdd_both_dACC_rse_gene <- means_mdd_both_dACC_rse_gene[which(names(means_mdd_both_dACC_rse_gene) %in% plasticity)] 																								
means_cont_both_dACC_rse_gene <- means_cont_both_dACC_rse_gene[which(names(means_cont_both_dACC_rse_gene) %in% plasticity)] 

means_ptsd_both_DLPFC_rse_gene <- means_ptsd_both_DLPFC_rse_gene[which(names(means_ptsd_both_DLPFC_rse_gene) %in% plasticity)] 
means_mdd_both_DLPFC_rse_gene <- means_mdd_both_DLPFC_rse_gene[which(names(means_mdd_both_DLPFC_rse_gene) %in% plasticity)] 																								
means_cont_both_DLPFC_rse_gene <- means_cont_both_DLPFC_rse_gene[which(names(means_cont_both_DLPFC_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dACC_rse_gene)

#seems meanExprs is rpkm (not subetting by disorder) and the above way is just counts (subsetting by disorder)
##meanExprs still not exactly mean of the disorders tho... 
#expression filter first?
##use log2(getRPKM) + 1?

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#Read in our data and assign unique name
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]


#Keep regions of interest
dACCkeepIndex = which(rse_gene$Region == "dACC")
dACC_rse_gene <- rse_gene[,dACCkeepIndex]
DLPFCkeepIndex = which(rse_gene$Region == "DLPFC")
DLPFC_rse_gene <- rse_gene[,DLPFCkeepIndex]

#Only keep samples where we have both the dACC and the DLPFC
test <- which(dACC_rse_gene$BrNum %in% DLPFC_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(DLPFC_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_DLPFC_rse_gene <- DLPFC_rse_gene[,test2]

#Check that the brain numbers are the same
both_DLPFC_rse_gene$BrNum <- as.character(both_DLPFC_rse_gene$BrNum)
unique(both_DLPFC_rse_gene$BrNum)
unique(both_dACC_rse_gene$BrNum)
identical(sort(both_DLPFC_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_DLPFC_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_DLPFC_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_DLPFC_rse_gene$Group == "PTSD")
ptsd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_DLPFC_rse_gene$Group == "MDD")
mdd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_DLPFC_rse_gene$Group == "Control")
cont_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_DLPFC_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_DLPFC_rse_gene)
dim(mdd_both_DLPFC_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(getRPKM(ptsd_both_dACC_rse_gene, "Length"))
means_mdd_both_dACC_rse_gene <- rowMeans(getRPKM(mdd_both_dACC_rse_gene, "Length"))
means_cont_both_dACC_rse_gene <- rowMeans(getRPKM(cont_both_dACC_rse_gene, "Length"))
means_cont_both_DLPFC_rse_gene <- rowMeans(getRPKM(cont_both_DLPFC_rse_gene, "Length"))
means_mdd_both_DLPFC_rse_gene <- rowMeans(getRPKM(mdd_both_DLPFC_rse_gene, "Length"))
means_ptsd_both_DLPFC_rse_gene <- rowMeans(getRPKM(ptsd_both_DLPFC_rse_gene, "Length"))

#Calculate Pearson correlations
cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dA
	
#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_DLPFC_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_DLPFC_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_DLPFC_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

#subset to genes of interest, then check correlation? 
#maybe look at exons, other features?

#ideas for genes to subset to
##BDNF, FKBP5, Mtor, Akt, Egr, NMDARs, AMPARs, CREB, NTRK2
###USF, Sp3/4
####mGluR2, 5HT2AR

#gencode IDs
plasticity <- c("ENSG00000176697.18","ENSG00000096060.14","ENSG00000198793.12","ENSG00000120738.7","ENSG00000164082.14","ENSG00000102468.10","ENSG00000176884.14","ENSG00000148053.15","ENSG00000118260.14")

means_ptsd_both_dACC_rse_gene <- means_ptsd_both_dACC_rse_gene[which(names(means_ptsd_both_dACC_rse_gene) %in% plasticity)] 
means_mdd_both_dACC_rse_gene <- means_mdd_both_dACC_rse_gene[which(names(means_mdd_both_dACC_rse_gene) %in% plasticity)] 																								
means_cont_both_dACC_rse_gene <- means_cont_both_dACC_rse_gene[which(names(means_cont_both_dACC_rse_gene) %in% plasticity)] 

means_ptsd_both_DLPFC_rse_gene <- means_ptsd_both_DLPFC_rse_gene[which(names(means_ptsd_both_DLPFC_rse_gene) %in% plasticity)] 
means_mdd_both_DLPFC_rse_gene <- means_mdd_both_DLPFC_rse_gene[which(names(means_mdd_both_DLPFC_rse_gene) %in% plasticity)] 																								
means_cont_both_DLPFC_rse_gene <- means_cont_both_DLPFC_rse_gene[which(names(means_cont_both_DLPFC_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dACC_rse_gene)


##use log2(getRPKM) + 1?

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#Read in our data and assign unique name
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]


#Keep regions of interest
dACCkeepIndex = which(rse_gene$Region == "dACC")
dACC_rse_gene <- rse_gene[,dACCkeepIndex]
DLPFCkeepIndex = which(rse_gene$Region == "DLPFC")
DLPFC_rse_gene <- rse_gene[,DLPFCkeepIndex]

#Only keep samples where we have both the dACC and the DLPFC
test <- which(dACC_rse_gene$BrNum %in% DLPFC_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(DLPFC_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_DLPFC_rse_gene <- DLPFC_rse_gene[,test2]

#Check that the brain numbers are the same
both_DLPFC_rse_gene$BrNum <- as.character(both_DLPFC_rse_gene$BrNum)
unique(both_DLPFC_rse_gene$BrNum)
unique(both_dACC_rse_gene$BrNum)
identical(sort(both_DLPFC_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_DLPFC_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_DLPFC_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_DLPFC_rse_gene$Group == "PTSD")
ptsd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_DLPFC_rse_gene$Group == "MDD")
mdd_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_DLPFC_rse_gene$Group == "Control")
cont_both_DLPFC_rse_gene <- both_DLPFC_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_DLPFC_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_DLPFC_rse_gene)
dim(mdd_both_DLPFC_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(log2(getRPKM(ptsd_both_dACC_rse_gene, "Length")+1))
means_mdd_both_dACC_rse_gene <- rowMeans(log2(getRPKM(mdd_both_dACC_rse_gene, "Length")+1))
means_cont_both_dACC_rse_gene <- rowMeans(log2(getRPKM(cont_both_dACC_rse_gene, "Length")+1))
means_cont_both_DLPFC_rse_gene <- rowMeans(log2(getRPKM(cont_both_DLPFC_rse_gene, "Length")+1))
means_mdd_both_DLPFC_rse_gene <- rowMeans(log2(getRPKM(mdd_both_DLPFC_rse_gene, "Length")+1))
means_ptsd_both_DLPFC_rse_gene <- rowMeans(log2(getRPKM(ptsd_both_DLPFC_rse_gene, "Length")+1))

#Calculate Pearson correlations
cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dACC_rse_gene)

#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_DLPFC_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_DLPFC_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_DLPFC_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

#subset to genes of interest, then check correlation? 
#maybe look at exons, other features?

#ideas for genes to subset to
##BDNF, FKBP5, Mtor, Akt, Egr, NMDARs, AMPARs, CREB, NTRK2
###USF, Sp3/4
####mGluR2, 5HT2AR

#gencode IDs
plasticity <- c("ENSG00000176697.18","ENSG00000096060.14","ENSG00000198793.12","ENSG00000120738.7","ENSG00000164082.14","ENSG00000102468.10","ENSG00000176884.14","ENSG00000148053.15","ENSG00000118260.14")

means_ptsd_both_dACC_rse_gene <- means_ptsd_both_dACC_rse_gene[which(names(means_ptsd_both_dACC_rse_gene) %in% plasticity)] 
means_mdd_both_dACC_rse_gene <- means_mdd_both_dACC_rse_gene[which(names(means_mdd_both_dACC_rse_gene) %in% plasticity)] 																								
means_cont_both_dACC_rse_gene <- means_cont_both_dACC_rse_gene[which(names(means_cont_both_dACC_rse_gene) %in% plasticity)] 

means_ptsd_both_DLPFC_rse_gene <- means_ptsd_both_DLPFC_rse_gene[which(names(means_ptsd_both_DLPFC_rse_gene) %in% plasticity)] 
means_mdd_both_DLPFC_rse_gene <- means_mdd_both_DLPFC_rse_gene[which(names(means_mdd_both_DLPFC_rse_gene) %in% plasticity)] 																								
means_cont_both_DLPFC_rse_gene <- means_cont_both_DLPFC_rse_gene[which(names(means_cont_both_DLPFC_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_DLPFC_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_DLPFC_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_DLPFC_rse_gene,means_cont_both_dACC_rse_gene)

#instead of taking means, leave as matrix?
cor(assay(ptsd_both_DLPFC_rse_gene), assay(ptsd_both_dACC_rse_gene))
#ensure cor is doing what i think its doing for all of the above analyses...

#maybe see how synaptic plasticity genes and the genes they regulate/affect are correlated across regions?



