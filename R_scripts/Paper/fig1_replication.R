###
## read back in merged stats
library(readxl)
library(SummarizedExperiment)

##########
## compare people
yale_pheno = read_excel("../../../LIBD Validation Cohort.xlsx")
load("../../../count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")

table(yale_pheno$BrNum %in% rse_gene$BrNum)

## our DEGs
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

## yale dacc only
mdd_yale = read.csv("../../../Yale_SuppTables_cleaned/S9-SupplementaryTable-MDD-DEGs-fixedgenenames.csv",
	as.is=TRUE,row.names=1)
mdd_yale_match = mdd_yale[match(deg_stats$ensemblID, mdd_yale$Geneid),]

ptsd_yale = read.csv("../../../Yale_SuppTables_cleaned/S1-SupplementaryTable-PTSD-DEGs-fixedgenenames.csv",
	as.is=TRUE,row.names=1)
ptsd_yale_match = ptsd_yale[match(deg_stats$ensemblID, ptsd_yale$Geneid),]

## yale full
ptsd_yale_full = read.csv("../../../Yale_SuppTables/EMAIL-Jaffe-DACC-PTSD-DEGs.csv",
	as.is=TRUE,row.names=1)
table(ptsd_yale$Geneid %in% ptsd_yale_full$Geneid)

ptsd_yale_full_match = ptsd_yale_full[
	match(deg_stats$ensemblID, ptsd_yale_full$Geneid),]

##########################
#### PTSD - LIBD signif ##
##########################

dacc_ptsd_sigIndex= which(deg_stats$dACC_adjPVal_PTSD < 0.1)
table(is.na(ptsd_yale_match$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex]))
# FALSE  TRUE
   # 57    15
table(is.na(ptsd_yale_full_match$log2FoldChange[dacc_ptsd_sigIndex]))
# FALSE  TRUE
   # 57    15

pdf("figures/figure1/libd_ptsd_dacc_replication_assessment.pdf")
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex],
	ptsd_yale_match$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex],
	pch=21,bg="grey",xlim = c(-1,1),ylim = c(-1,1), 
	main = "PTSD vs Control: log2 FCs",
	xlab = "Discovery Dataset", ylab = "Replication Dataset")
abline(v=0,h=0, lty=2)
cor.test(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex],
	ptsd_yale_match$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex])
dev.off()

# t = 5.3822, df = 55, p-value = 1.564e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.3858456 0.7353804
# sample estimates:
      # cor
# 0.5873607

table(sign(ptsd_yale_match$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex]))
# FALSE  TRUE
   # 11    46

table(sign(ptsd_yale_match$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex]) & 
  ptsd_yale_match$PTSD.dACC.padj[dacc_ptsd_sigIndex] < 0.1)
# FALSE  TRUE
   # 48     9

## full table
table(sign(ptsd_yale_full_match$log2FoldChange[dacc_ptsd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex]) & 
  ptsd_yale_full_match$pvalue[dacc_ptsd_sigIndex] < 0.05)
 # FALSE  TRUE
   # 45    12

table(sign(ptsd_yale_full_match$log2FoldChange[dacc_ptsd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_PTSD[dacc_ptsd_sigIndex]) & 
  ptsd_yale_full_match$pvalue[dacc_ptsd_sigIndex] < 0.1)
 # FALSE  TRUE
#	   37    20

##########################
#### MDD - LIBD signif ##
##########################

dacc_mdd_sigIndex= which(deg_stats$dACC_adjPVal_MDD < 0.1)
table(is.na(mdd_yale_match$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex]))
# FALSE  TRUE
  # 202    47

pdf("figures/figure1/libd_mdd_dacc_replication_assessment.pdf")
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(deg_stats$dACC_logFC_MDD[dacc_mdd_sigIndex],
	mdd_yale_match$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex],
	pch=21,bg="grey",xlim = c(-1,1),ylim = c(-1,1), 
	main = "MDD vs Control: log2 FCs",
	xlab = "Discovery Dataset", ylab = "Replication Dataset")
cor.test(deg_stats$dACC_logFC_MDD[dacc_mdd_sigIndex],
	mdd_yale_match$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex])
	abline(v=0,h=0, lty=2)

dev.off()
# t = 7.4508, df = 200, p-value = 2.731e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.3506310 0.5676436
# sample estimates:
      # cor
# 0.4661195


table(sign(mdd_yale_match$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_MDD[dacc_mdd_sigIndex]))
# FALSE  TRUE
   # 64   138

table(sign(mdd_yale_match$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex]) ==
 sign(deg_stats$dACC_logFC_MDD[dacc_mdd_sigIndex]) & 
  mdd_yale_match$MDD.dACC.padj[dacc_mdd_sigIndex] < 0.1)
# FALSE  TRUE
  # 195     7


##########################
#### PTSD - Yale signif ##
##########################

dacc_ptsd_sigIndex_yale = which(ptsd_yale$PTSD.dACC.padj < 0.1)
deg_stats_match = deg_stats[match(ptsd_yale$Geneid, deg_stats$ensemblID),]
table(is.na(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]))
# FALSE  TRUE
  # 193     3


plot(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale],
	ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale])
cor.test(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale],
	ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale])

# t = 9.3613, df = 191, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.4556668 0.6505262
# sample estimates:
      # cor
# 0.5608138


table(sign(ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]))
# FALSE  TRUE
   # 34   159

table(sign(ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_PTSD[dacc_ptsd_sigIndex_yale] < 0.1)
# FALSE  TRUE
  # 112    81

table(sign(ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_PTSD[dacc_ptsd_sigIndex_yale] < 0.05)
# FALSE  TRUE
#   136    57
table(sign(ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_PTSD[dacc_ptsd_sigIndex_yale] < 0.01)
# FALSE  TRUE
#   166    27
table(sign(ptsd_yale$PTSD.dACC.log2FoldChange[dacc_ptsd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_PTSD[dacc_ptsd_sigIndex_yale]) & 
  deg_stats_match$dACC_adjPVal_PTSD[dacc_ptsd_sigIndex_yale] < 0.1)
# FALSE  TRUE
#   184     9

##########################
#### MDD - Yale signif ##
##########################

dacc_mdd_sigIndex_yale = which(mdd_yale$MDD.dACC.padj < 0.1)
deg_stats_match = deg_stats[match(mdd_yale$Geneid, deg_stats$ensemblID),]
table(is.na(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]))
# FALSE  TRUE
#   182     7


plot(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale],
	mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale])
cor.test(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale],
	mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale])

# t = 6.1995, df = 180, p-value = 3.768e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.2918177 0.5324377
# sample estimates:
      # cor
# 0.4194683



table(sign(mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]))
# FALSE  TRUE
#    59   123

table(sign(mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_MDD[dacc_mdd_sigIndex_yale] < 0.1)
  
table(sign(mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_MDD[dacc_mdd_sigIndex_yale] < 0.05)
# FALSE  TRUE
#   155    27
table(sign(mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]) & 
  deg_stats_match$dACC_PValue_MDD[dacc_mdd_sigIndex_yale] < 0.01)
# FALSE  TRUE
#   169    13
table(sign(mdd_yale$MDD.dACC.log2FoldChange[dacc_mdd_sigIndex_yale]) ==
 sign(deg_stats_match$dACC_logFC_MDD[dacc_mdd_sigIndex_yale]) & 
  deg_stats_match$dACC_adjPVal_MDD[dacc_mdd_sigIndex_yale] < 0.1)
# FALSE  TRUE
#   175     7
