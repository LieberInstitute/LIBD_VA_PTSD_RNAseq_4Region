library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(recount)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load("rdas/mergedEqtl_output_4regionsmerged_allFeatures_final_withfdr.rda")
load("rdas/mergedEqtl_output_BasoAmyg_allFeatures_annotation.rda")
load("rdas/mergedEqtl_output_MedialAmyg_allFeatures_annotation.rda")
load("rdas/mergedEqtl_output_dACC_allFeatures_annotation.rda")
load("rdas/mergedEqtl_output_DLPFC_allFeatures_annotation.rda")


supp3 <- read.csv(file="csvs/41593_2019_447_MOESM3_ESM.csv", header=TRUE)
head(supp3)
class(supp3)
supp4 <- read.csv(file="csvs/41593_2019_447_MOESM4_ESM.csv", header=TRUE)
ls()

BLA_annotated$rsSNPID <- ss(BLA_annotated$ID, ":", 1)
MeA_annotated$rsSNPID <- ss(MeA_annotated$ID, ":", 1)
dACC_annotated$rsSNPID <- ss(dACC_annotated$ID, ":", 1)
DLPFC_annotated$rsSNPID <- ss(DLPFC_annotated$ID, ":", 1)


#do all the snps in BLA_annotated exist in supp3 and vice versa?
summary(BLA_annotated$rsSNPID %in% supp3$SNP)
#    Mode    FALSE     TRUE 
# logical 19704375    11956 
summary(supp3$SNP %in% BLA_annotated$rsSNPID )
#   Mode   FALSE    TRUE 
#logical    8858    1142 
summary(unique(BLA_annotated$rsSNPID) %in% supp3$SNP)
#   Mode   FALSE    TRUE 
#logical 3513200    1142 

#subset
BLA_subset <- subset(BLA_annotated, rsSNPID %in% supp3$SNP)