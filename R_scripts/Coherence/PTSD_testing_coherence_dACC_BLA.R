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
BasoAmygkeepIndex = which(rse_gene$Region == "BasoAmyg")
BasoAmyg_rse_gene <- rse_gene[,BasoAmygkeepIndex]

#Only keep samples where we have both the dACC and the BasoAmyg
test <- which(dACC_rse_gene$BrNum %in% BasoAmyg_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(BasoAmyg_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_BasoAmyg_rse_gene <- BasoAmyg_rse_gene[,test2]

#Check that the brain numbers are the same
both_BasoAmyg_rse_gene$BrNum <- as.character(both_BasoAmyg_rse_gene$BrNum)
unique(both_BasoAmyg_rse_gene$BrNum)
unique(both_dACC_rse_gene$BrNum)
identical(sort(both_BasoAmyg_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_BasoAmyg_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_BasoAmyg_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "PTSD")
ptsd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "MDD")
mdd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "Control")
cont_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_BasoAmyg_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_BasoAmyg_rse_gene)
dim(mdd_both_BasoAmyg_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(assay(ptsd_both_dACC_rse_gene))
means_mdd_both_dACC_rse_gene <- rowMeans(assay(mdd_both_dACC_rse_gene))
means_cont_both_dACC_rse_gene <- rowMeans(assay(cont_both_dACC_rse_gene))
means_cont_both_BasoAmyg_rse_gene <- rowMeans(assay(cont_both_BasoAmyg_rse_gene))
means_mdd_both_BasoAmyg_rse_gene <- rowMeans(assay(mdd_both_BasoAmyg_rse_gene))
means_ptsd_both_BasoAmyg_rse_gene <- rowMeans(assay(ptsd_both_BasoAmyg_rse_gene))

#Calculate Pearson correlations
cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dACC_rse_gene)

#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_BasoAmyg_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_BasoAmyg_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_BasoAmyg_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

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

means_ptsd_both_BasoAmyg_rse_gene <- means_ptsd_both_BasoAmyg_rse_gene[which(names(means_ptsd_both_BasoAmyg_rse_gene) %in% plasticity)] 
means_mdd_both_BasoAmyg_rse_gene <- means_mdd_both_BasoAmyg_rse_gene[which(names(means_mdd_both_BasoAmyg_rse_gene) %in% plasticity)] 																								
means_cont_both_BasoAmyg_rse_gene <- means_cont_both_BasoAmyg_rse_gene[which(names(means_cont_both_BasoAmyg_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dACC_rse_gene)

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
BasoAmygkeepIndex = which(rse_gene$Region == "BasoAmyg")
BasoAmyg_rse_gene <- rse_gene[,BasoAmygkeepIndex]

#Only keep samples where we have both the dACC and the BasoAmyg
test <- which(dACC_rse_gene$BrNum %in% BasoAmyg_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(BasoAmyg_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_BasoAmyg_rse_gene <- BasoAmyg_rse_gene[,test2]

#Check that the brain numbers are the same
both_BasoAmyg_rse_gene$BrNum <- as.character(both_BasoAmyg_rse_gene$BrNum)
unique(both_BasoAmyg_rse_gene$BrNum)
unique(both_dACC_rse_gene$BrNum)
identical(sort(both_BasoAmyg_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_BasoAmyg_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_BasoAmyg_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "PTSD")
ptsd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "MDD")
mdd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "Control")
cont_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_BasoAmyg_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_BasoAmyg_rse_gene)
dim(mdd_both_BasoAmyg_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(getRPKM(ptsd_both_dACC_rse_gene, "Length"))
means_mdd_both_dACC_rse_gene <- rowMeans(getRPKM(mdd_both_dACC_rse_gene, "Length"))
means_cont_both_dACC_rse_gene <- rowMeans(getRPKM(cont_both_dACC_rse_gene, "Length"))
means_cont_both_BasoAmyg_rse_gene <- rowMeans(getRPKM(cont_both_BasoAmyg_rse_gene, "Length"))
means_mdd_both_BasoAmyg_rse_gene <- rowMeans(getRPKM(mdd_both_BasoAmyg_rse_gene, "Length"))
means_ptsd_both_BasoAmyg_rse_gene <- rowMeans(getRPKM(ptsd_both_BasoAmyg_rse_gene, "Length"))

#Calculate Pearson correlations
cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dA
	
#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_BasoAmyg_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_BasoAmyg_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_BasoAmyg_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

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

means_ptsd_both_BasoAmyg_rse_gene <- means_ptsd_both_BasoAmyg_rse_gene[which(names(means_ptsd_both_BasoAmyg_rse_gene) %in% plasticity)] 
means_mdd_both_BasoAmyg_rse_gene <- means_mdd_both_BasoAmyg_rse_gene[which(names(means_mdd_both_BasoAmyg_rse_gene) %in% plasticity)] 																								
means_cont_both_BasoAmyg_rse_gene <- means_cont_both_BasoAmyg_rse_gene[which(names(means_cont_both_BasoAmyg_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dACC_rse_gene)


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
BasoAmygkeepIndex = which(rse_gene$Region == "BasoAmyg")
BasoAmyg_rse_gene <- rse_gene[,BasoAmygkeepIndex]

#Only keep samples where we have both the dACC and the BasoAmyg
test <- which(dACC_rse_gene$BrNum %in% BasoAmyg_rse_gene$BrNum)
both_dACC_rse_gene <- dACC_rse_gene[,test]
test2 <- which(BasoAmyg_rse_gene$BrNum %in% both_dACC_rse_gene$BrNum)
both_BasoAmyg_rse_gene <- BasoAmyg_rse_gene[,test2]

#Check that the brain numbers are the same
both_BasoAmyg_rse_gene$BrNum <- as.character(both_BasoAmyg_rse_gene$BrNum)
length(unique(both_BasoAmyg_rse_gene$BrNum))
length(unique(both_dACC_rse_gene$BrNum))
identical(sort(both_BasoAmyg_rse_gene$BrNum), sort(both_dACC_rse_gene$BrNum))

identical(unique(both_BasoAmyg_rse_gene$BrNum), unique(both_dACC_rse_gene$BrNum))


#Subset into the different diagnoses
summary(as.factor(both_dACC_rse_gene$Group))
summary(as.factor(both_BasoAmyg_rse_gene$Group))

ptsdkeepIndex <- which(both_dACC_rse_gene$Group == "PTSD")
ptsd_both_dACC_rse_gene <- both_dACC_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_dACC_rse_gene$Group == "MDD")
mdd_both_dACC_rse_gene <- both_dACC_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_dACC_rse_gene$Group == "Control")
cont_both_dACC_rse_gene <- both_dACC_rse_gene[,contkeepIndex]

ptsdkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "PTSD")
ptsd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,ptsdkeepIndex]
MDDkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "MDD")
mdd_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,MDDkeepIndex]
contkeepIndex <- which(both_BasoAmyg_rse_gene$Group == "Control")
cont_both_BasoAmyg_rse_gene <- both_BasoAmyg_rse_gene[,contkeepIndex]

#Check that the number of samples in each diagnosis match
dim(cont_both_BasoAmyg_rse_gene)
dim(cont_both_dACC_rse_gene)
dim(ptsd_both_dACC_rse_gene)
dim(ptsd_both_BasoAmyg_rse_gene)
dim(mdd_both_BasoAmyg_rse_gene)
dim(mdd_both_dACC_rse_gene)

#Get the average expression for each gene across samples in each diagnosis
means_ptsd_both_dACC_rse_gene <- rowMeans(log2(getRPKM(ptsd_both_dACC_rse_gene, "Length")+1))
means_mdd_both_dACC_rse_gene <- rowMeans(log2(getRPKM(mdd_both_dACC_rse_gene, "Length")+1))
means_cont_both_dACC_rse_gene <- rowMeans(log2(getRPKM(cont_both_dACC_rse_gene, "Length")+1))
means_cont_both_BasoAmyg_rse_gene <- rowMeans(log2(getRPKM(cont_both_BasoAmyg_rse_gene, "Length")+1))
means_mdd_both_BasoAmyg_rse_gene <- rowMeans(log2(getRPKM(mdd_both_BasoAmyg_rse_gene, "Length")+1))
means_ptsd_both_BasoAmyg_rse_gene <- rowMeans(log2(getRPKM(ptsd_both_BasoAmyg_rse_gene, "Length")+1))

#Calculate Pearson correlations
cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dACC_rse_gene)

#next steps

#only evaluate genes that are in both?
summary(names(means_ptsd_both_BasoAmyg_rse_gene) %in% names(means_ptsd_both_dACC_rse_gene))
summary(names(means_mdd_both_BasoAmyg_rse_gene) %in% names(means_mdd_both_dACC_rse_gene))
summary(names(means_cont_both_BasoAmyg_rse_gene) %in% names(means_cont_both_dACC_rse_gene))

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

means_ptsd_both_BasoAmyg_rse_gene <- means_ptsd_both_BasoAmyg_rse_gene[which(names(means_ptsd_both_BasoAmyg_rse_gene) %in% plasticity)] 
means_mdd_both_BasoAmyg_rse_gene <- means_mdd_both_BasoAmyg_rse_gene[which(names(means_mdd_both_BasoAmyg_rse_gene) %in% plasticity)] 																								
means_cont_both_BasoAmyg_rse_gene <- means_cont_both_BasoAmyg_rse_gene[which(names(means_cont_both_BasoAmyg_rse_gene) %in% plasticity)] 

cor(means_ptsd_both_BasoAmyg_rse_gene,means_ptsd_both_dACC_rse_gene)
cor(means_mdd_both_BasoAmyg_rse_gene,means_mdd_both_dACC_rse_gene)
cor(means_cont_both_BasoAmyg_rse_gene,means_cont_both_dACC_rse_gene)

#instead of taking means, leave as matrix?
cor(assay(ptsd_both_BasoAmyg_rse_gene), assay(ptsd_both_dACC_rse_gene))
#ensure cor is doing what i think its doing for all of the above analyses...

#maybe see how synaptic plasticity genes and the genes they regulate/affect are correlated across regions?


###########################################################################################################################
##use log2(getRPKM) + 0.5 per individual (compute individual level correlations across features) 
#based on github: `brainseq_phase2/correlation/corr_analysis.R` and ``brainseq_phase2/correlation/corr_individual.R`
###########################################################################################################################

#Read in our data
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')
library(sva)
library('readxl')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])
#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])
#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])
#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])

## expression filter to remove lowly expressed stuff (note, need like 200Gb of RAM for jxn filtering)
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]
#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]
#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]
#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]

#load qSVA model, then save with keepIndex for each of the regions of interest. This will help with "apply" functions later.
#Note, the qSVs were geneerated across all four regions (dACC, DLPFC, BLA, and MeA) at once with "region" in the model (see `qSV_model_anaysis_BasoAmyg_PTSD_updated.R`). 
#The "keepIndex" is the only part specific to a region.

load('rdas/PTSD_qsvs.Rdata')
keepIndex = which(rse_gene$Region == "BasoAmyg")
save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'rdas/PTSD_qsvs_BasoAmyg.Rdata')	
	
load('rdas/PTSD_qsvs.Rdata')
keepIndex = which(rse_gene$Region == "dACC")
save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'rdas/PTSD_qsvs_dACC.Rdata')

regions <- c('BasoAmyg', 'dACC')
qinfo <- lapply(regions, function(region) {
    f <- paste0('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/Coherence/PTSD_qsvs_', region, '.Rdata')
    message(paste(Sys.time(), 'loading file', f))
    load(f, verbose = TRUE)
    res <- list(
        'qsvBonf' = qsvBonf,
        'qSVs' = qSVs,
        'mod' = mod,
        'modQsva' = modQsva,
        'keepIndex' = keepIndex
    )
    return(res)
})
names(qinfo) <- regions

## Keep only samples that were observed twice
brains <- unlist(lapply(regions, function(region) {
    colData(rse_gene)$BrNum[qinfo[[region]]$keepIndex]
}))
both <- brains[which(duplicated(brains))]
length(both)
#320

## Subset the rses for each brain region and re-order
simple_rse <- lapply(regions, function(region) {
    
    res <- list(
        'gene' = rse_gene[, qinfo[[region]]$keepIndex],
        'exon' = rse_exon[, qinfo[[region]]$keepIndex],
        'jxn' = rse_jx[, qinfo[[region]]$keepIndex],
        'tx' = rse_tx[, qinfo[[region]]$keepIndex]
    )
    
    m <- match(both, colData(res$gene)$BrNum)
    stopifnot(all(!is.na(m)))
    ## Re-order by brain id
    lapply(res, function(x) {
        x[, m]
    })
})
names(simple_rse) <- regions

## Keep only dup brains on the mods and re-order
modQsva <- lapply(regions, function(region) {
    m <- match(both, colData(rse_gene[, qinfo[[region]]$keepIndex])$BrNum)
    stopifnot(all(!is.na(m)))
    
    qinfo[[region]]$modQsva[m, ]
})
names(modQsva) <- regions
stopifnot(all(sapply(modQsva, nrow) == 320))
stopifnot(all(sapply(modQsva, function(x) { colnames(x)[2] }) == 'GroupMDD'))
stopifnot(all(sapply(modQsva, function(x) { colnames(x)[3] }) == 'GroupPTSD'))
stopifnot(all(sapply(modQsva, function(x) { colnames(x)[4] }) == 'AgeDeath'))
stopifnot(all(sapply(modQsva, function(x) { colnames(x)[5] }) == 'SexM'))

message(paste(Sys.time(), 'saving rse and modQsva info'))
save(simple_rse, modQsva, file = 'rdas/Coherence/rse_and_modQsva_BasoAmyg_dACC.Rdata')

## Extract expr
#note, we havent subset the rse_gene object to the disorder of interest so this isnt specific to any disorder
expr <- lapply(simple_rse, function(rses) {
    res <- list(
        "geneRpkm" = getRPKM(rses$gene,"Length"),
        "exonRpkm" = getRPKM(rses$exon,"Length"),
        "jxnRp10m" = getRPKM(rses$jxn, "Length"),
        "txTpm" = assays(rses$tx)$tpm
    )
    lapply(res, function(x) { log2(x + 0.5) })
})

save(expr, file = 'rdas/Coherence/expr_BasoAmyg_dACC.Rdata')

## Create cleaned expr versions protecting Dx, age, and sex

#time to subset by disorder
#here, we protect the MDD diagnosis variable below so this cleaned expression is for MDD
modQsva <- lapply(regions,function(region) {
	modQsva[[region]][,c(1,2,4,5,3,6:28)]
})			   
names(modQsva) <- regions			 

cleaned_MDD <- lapply(regions, function(region) {
    lapply(expr[[region]], cleaningY, modQsva[[region]], P = 4)
})
names(cleaned_MDD) <- regions
sapply(cleaned_MDD, function(x) sapply(x, dim) )
save(cleaned_MDD, file = 'rdas/Coherence/cleaned_BasoAmyg_dACC_MDD.Rdata')

#subset the simple_rse, expr object, and cleaned object for MDD and controls only (makes plotting easier later)
#note, just collecting index from BasoAmyg because both regions and all features of each region are already in same order (if they werent, you couldnt do it this way)	   

#simple rse
types <- names(simple_rse[['BasoAmyg']])
Cont_MDDIndex <- which(colData(simple_rse[['BasoAmyg']][['gene']])$Group == "MDD" | colData(simple_rse[['BasoAmyg']][['gene']])$Group == "Control")

simple_rse_ContMDD <- lapply(regions, function(region) {
		lapply(types, function(type) {
		simple_rse[[region]][[type]][,Cont_MDDIndex]
			})
})
names(simple_rse_ContMDD) <- regions
names(simple_rse_ContMDD[[1]]) <-  types
names(simple_rse_ContMDD[[2]]) <-  types

#expr object
types2 <- names(expr[[1]])
expr_ContMDD <- lapply(regions, function(region) {
		lapply(types2, function(type) {
		expr[[region]][[type]][,Cont_MDDIndex]
			})
})
names(expr_ContMDD) <- regions
names(expr_ContMDD[[1]]) <-  types2
names(expr_ContMDD[[2]]) <-  types2

#cleaned expr object
types3 <- names(cleaned_MDD[[1]])
cleaned_ContMDD <- lapply(regions, function(region) {
		lapply(types3, function(type) {
		cleaned_MDD[[region]][[type]][,Cont_MDDIndex]
			})
})
names(cleaned_ContMDD) <- regions
names(cleaned_ContMDD[[1]]) <-  types3
names(cleaned_ContMDD[[2]]) <-  types3

save(expr_ContMDD, file = 'rdas/Coherence/expr_BasoAmyg_dACC_MDD_subset.Rdata')
save(cleaned_ContMDD, file = 'rdas/Coherence/cleaned_BasoAmyg_dACC_MDD_subset.Rdata')

## Function for computing paired cors
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(ncol(x)), function(i) {
        cor(x[, i], y[, i])
    })
}

## Test paircor
test <- cor(cleaned_ContMDD[['BasoAmyg']][['geneRpkm']][1:100, ], cleaned_ContMDD[['dACC']][['geneRpkm']][1:100, ])
#test is subject vs subject matrix, 320x320
test2 <- paircor(cleaned_ContMDD[['BasoAmyg']][['geneRpkm']][1:100, ], cleaned_ContMDD[['dACC']][['geneRpkm']][1:100, ])
#test2 has length 320 and is correlation values of correlation between the two brain regions for each individual
stopifnot(all(diag(test) - test2 == 0))

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['BasoAmyg']][[feature]], exp[['dACC']][[feature]])
    })
    names(res) <- sets
    return(res)
}

## Compute individual level correlations across features
indv_expr_ContMDD <- computecor(expr_ContMDD)
indv_cleaned_ContMDD <- computecor(cleaned_ContMDD)
save(indv_expr_ContMDD, indv_cleaned_ContMDD, file = 'rdas/Coherence/indv_corr_BasoAmyg_dACC_MDD_subset.Rdata')


#PTSD next...
#Read in our data
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')
library(sva)
library('readxl')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')
load('rdas/Coherence/rse_and_modQsva_BasoAmyg_dACC.Rdata',verbose=TRUE)
load('rdas/Coherence/expr_BasoAmyg_dACC.Rdata',verbose=TRUE)

regions <- c('BasoAmyg', 'dACC')

## Create cleaned expr versions protecting Dx, age, and sex
#here we subset by disorder
#here, we protect the PTSD diagnosis variable below so this cleaned expression is for PTSD
modQsva <- lapply(regions,function(region) {
	modQsva[[region]][,c(1,3,4,5,2,6:28)]
})			   
names(modQsva) <- regions			 

cleaned_PTSD <- lapply(regions, function(region) {
    lapply(expr[[region]], cleaningY, modQsva[[region]], P = 4)
})
names(cleaned_PTSD) <- regions
sapply(cleaned_PTSD, function(x) sapply(x, dim) )
save(cleaned_PTSD, file = 'rdas/Coherence/cleaned_BasoAmyg_dACC_PTSD.Rdata')

#subset the simple_rse, expr object, and cleaned object for PTSD and controls only (makes plotting easier later)
#note, just collecting index from BasoAmyg because both regions and all features of each region are already in same order (if they werent, you couldnt do it this way)	   

#simple rse
types <- names(simple_rse[['BasoAmyg']])
Cont_PTSDIndex <- which(colData(simple_rse[['BasoAmyg']][['gene']])$Group == "PTSD" | colData(simple_rse[['BasoAmyg']][['gene']])$Group == "Control")

simple_rse_ContPTSD <- lapply(regions, function(region) {
		lapply(types, function(type) {
		simple_rse[[region]][[type]][,Cont_PTSDIndex]
			})
})
names(simple_rse_ContPTSD) <- regions
names(simple_rse_ContPTSD[[1]]) <-  types
names(simple_rse_ContPTSD[[2]]) <-  types

#expr object
types2 <- names(expr[[1]])
expr_ContPTSD <- lapply(regions, function(region) {
		lapply(types2, function(type) {
		expr[[region]][[type]][,Cont_PTSDIndex]
			})
})
names(expr_ContPTSD) <- regions
names(expr_ContPTSD[[1]]) <-  types2
names(expr_ContPTSD[[2]]) <-  types2

#cleaned expr object
types3 <- names(cleaned_PTSD[[1]])
cleaned_ContPTSD <- lapply(regions, function(region) {
		lapply(types3, function(type) {
		cleaned_PTSD[[region]][[type]][,Cont_PTSDIndex]
			})
})
names(cleaned_ContPTSD) <- regions
names(cleaned_ContPTSD[[1]]) <-  types3
names(cleaned_ContPTSD[[2]]) <-  types3

save(expr_ContPTSD, file = 'rdas/Coherence/expr_BasoAmyg_dACC_PTSD_subset.Rdata')
save(cleaned_ContPTSD, file = 'rdas/Coherence/cleaned_BasoAmyg_dACC_PTSD_subset.Rdata')

## Function for computing paired cors
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(ncol(x)), function(i) {
        cor(x[, i], y[, i])
    })
}

## Test paircor
test <- cor(cleaned_ContPTSD[['BasoAmyg']][['geneRpkm']][1:100, ], cleaned_ContPTSD[['dACC']][['geneRpkm']][1:100, ])
#test is subject vs subject matrix
test2 <- paircor(cleaned_ContPTSD[['BasoAmyg']][['geneRpkm']][1:100, ], cleaned_ContPTSD[['dACC']][['geneRpkm']][1:100, ])
#test2 is correlation values of correlation between the two brain regions for each individual
stopifnot(all(diag(test) - test2 == 0))

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['BasoAmyg']][[feature]], exp[['dACC']][[feature]])
    })
    names(res) <- sets
    return(res)
}

## Compute individual level correlations across features
indv_expr_ContPTSD <- computecor(expr_ContPTSD)
indv_cleaned_ContPTSD <- computecor(cleaned_ContPTSD)
save(indv_expr_ContPTSD, indv_cleaned_ContPTSD, file = 'rdas/Coherence/indv_corr_BasoAmyg_dACC_PTSD_subset.Rdata')
	   

################################
#Plotting of individual correlations

#plot MDD
pdf('pdf/Coherence/pdf/indv_box_BasoAmyg_dACC_MDD.pdf', useDingbats = FALSE)
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Group[m]
    dx <- factor(ifelse(dx == 'MDD', 'MDD', ifelse(dx == 'Control', 'Control', 'hmm')))
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4] #this should give p value
    ylim <- range(x[m])
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'MDD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = 'light blue')
}, indv_expr_ContMDD, simple_rse_ContMDD[[1]], names(indv_expr_ContMDD), 'expr')
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Group[m]
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    dx <- factor(ifelse(dx == 'MDD', 'MDD', ifelse(dx == 'Control', 'Control', 'hmm')))
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'MDD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = '#009E73')
}, indv_cleaned_ContMDD, simple_rse_ContMDD[[1]], names(indv_cleaned_ContMDD), 'cleaned expr (keeping Dx, Sex, Age)')
dev.off()


#plot PTSD
pdf('pdf/Coherence/indv_box_BasoAmyg_dACC_PTSD.pdf', useDingbats = FALSE)
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Group[m]
    dx <- factor(ifelse(dx == 'PTSD', 'PTSD', ifelse(dx == 'Control', 'Control', 'hmm')))
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4] #this should give p value
    ylim <- range(x[m])
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'PTSD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = 'light blue')
}, indv_expr_ContPTSD, simple_rse_ContPTSD[[1]], names(indv_expr_ContPTSD), 'expr')
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Group[m]
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    dx <- factor(ifelse(dx == 'PTSD', 'PTSD', ifelse(dx == 'Control', 'Control', 'hmm')))
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'PTSD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = '#009E73')
}, indv_cleaned_ContPTSD, simple_rse_ContPTSD[[1]], names(indv_cleaned_ContPTSD), 'cleaned expr (keeping Dx, Sex, Age)')
dev.off()



## Explore covariates to see if low correlation is related to any main covariate
#note, our cleaned objects take into account age and sex
#note the origina DEG model was:
#mod = model.matrix(~Group + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate 
#	+ ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

#MDD
pdf('pdf/Coherence/indv_corr_vs_covariates_BasoAmyg_dACC_MDD.pdf', useDingbats = FALSE)
lapply(c('AgeDeath'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse_ContMDD[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_expr_ContMDD, names(indv_expr_ContMDD))
		mapply(function(x, set) {
           	plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_cleaned_ContMDD, names(indv_cleaned_ContMDD))
        return(NULL)
})
lapply(c('Sex'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse_ContMDD[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr')
        }, indv_expr_ContMDD, names(indv_expr_ContMDD))
		mapply(function(x, set) {
           	plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned')
        }, indv_cleaned_ContMDD, names(indv_cleaned_ContMDD))
        return(NULL)
})

lapply(c('snpPC1', 'snpPC2', 'snpPC3', 'snpPC8', 'snpPC9', 'snpPC10', 'RIN', 'rRNA_rate', 'totalAssignedGene', 'mitoRate', 'ERCCsumLogErr', 'overallMapRate'), function(var) {
    message(paste(Sys.time(), 'processing', var))
    lapply(names(simple_rse_ContMDD), function(region) {
        pd <- colData(simple_rse_ContMDD[[region]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - expr', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_expr_ContMDD, names(indv_expr_ContMDD))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - cleaned', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_cleaned_ContMDD, names(indv_cleaned_ContMDD))
        return(NULL)
    })
})
dev.off()

#PTSD
pdf('pdf/Coherence/indv_corr_vs_covariates_BasoAmyg_dACC_PTSD.pdf', useDingbats = FALSE)
lapply(c('AgeDeath'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse_ContPTSD[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_expr_ContPTSD, names(indv_expr_ContPTSD))
		mapply(function(x, set) {
           	plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_cleaned_ContPTSD, names(indv_cleaned_ContPTSD))
        return(NULL)
})
lapply(c('Sex'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse_ContPTSD[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr')
        }, indv_expr_ContPTSD, names(indv_expr_ContPTSD))
		mapply(function(x, set) {
           	plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned')
        }, indv_cleaned_ContPTSD, names(indv_cleaned_ContPTSD))
        return(NULL)
})

lapply(c('snpPC1', 'snpPC2', 'snpPC3', 'snpPC8', 'snpPC9', 'snpPC10', 'RIN', 'rRNA_rate', 'totalAssignedGene', 'mitoRate', 'ERCCsumLogErr', 'overallMapRate'), function(var) {
    message(paste(Sys.time(), 'processing', var))
    lapply(names(simple_rse_ContPTSD), function(region) {
        pd <- colData(simple_rse_ContPTSD[[region]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_expr_ContPTSD, names(indv_expr_ContPTSD))
		mapply(function(x, set) {
           	plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned', pch = 19, col = as.factor(pd$Group))
			legend('topleft',unique(pd$Group),col=1:length(unique(pd$Group)),pch=19)
        }, indv_cleaned_ContPTSD, names(indv_cleaned_ContPTSD))
        return(NULL)
    })
})
dev.off()
