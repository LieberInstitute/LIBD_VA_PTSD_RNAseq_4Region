library(SummarizedExperiment)
library(readxl)
library(jaffelab)
library(ggplot2)
library(gplots)
library(RColorBrewer)


setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')


#MDD
li_BasoAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_MDD_Li2018.csv")
li_MedialAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_MDD_Li2018.csv")
li_dACC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_Li2018.csv")
li_DLPFC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_MDD_Li2018.csv")


lake_BasoAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_MDD_Lake2018.csv")
lake_MedialAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_MDD_Lake2018.csv")
lake_dACC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_Lake2018.csv")
lake_DLPFC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_MDD_Lake2018.csv")


lake_broad_BasoAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_MDD_Lakebroad2018.csv")
lake_broad_MedialAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_MDD_Lakebroad2018.csv")
lake_broad_dACC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_Lakebroad2018.csv")
lake_broad_DLPFC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_MDD_Lakebroad2018.csv")

colnames(li_DLPFC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_dACC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_MedialAmyg) <-c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_BasoAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")

colnames(lake_DLPFC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_dACC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_MedialAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_BasoAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")


colnames(lake_broad_DLPFC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_broad_dACC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_broad_MedialAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(lake_broad_BasoAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")


#Li
rnames <-li_BasoAmyg[,1]
mat <- data.matrix(li_BasoAmyg[,2:ncol(li_BasoAmyg)])
rownames(mat) <- rnames

pdf("pdf/BLA_CellType_MDD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "BLA MDD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_MedialAmyg[,1]
mat <- data.matrix(li_MedialAmyg[,2:ncol(li_MedialAmyg)])
rownames(mat) <- rnames

pdf("pdf/MedialAmyg_CellType_MDD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "MeA MDD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_DLPFC[,1]
mat <- data.matrix(li_DLPFC[,2:ncol(li_DLPFC)])
rownames(mat) <- rnames

pdf("pdf/DLPFC_CellType_MDD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "DLPFC MDD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_dACC[,1]
mat <- data.matrix(li_dACC[,2:ncol(li_dACC)])
rownames(mat) <- rnames

pdf("pdf/dACC_CellType_MDD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "dACC MDD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device




###############################################################################

li_BasoAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_PTSD_Li2018.csv")
li_MedialAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_PTSD_Li2018.csv")
li_dACC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_Li2018.csv")
li_DLPFC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_PTSD_Li2018.csv")

colnames(li_DLPFC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_dACC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_MedialAmyg) <-c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_BasoAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")

#PTSD
#Li
rnames <-li_BasoAmyg[,1]
mat <- data.matrix(li_BasoAmyg[,2:ncol(li_BasoAmyg)])
rownames(mat) <- rnames

pdf("pdf/BLA_CellType_PTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "BLA PTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_MedialAmyg[,1]
mat <- data.matrix(li_MedialAmyg[,2:ncol(li_MedialAmyg)])
rownames(mat) <- rnames

pdf("pdf/MedialAmyg_CellType_PTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "MeA PTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_DLPFC[,1]
mat <- data.matrix(li_DLPFC[,2:ncol(li_DLPFC)])
rownames(mat) <- rnames

pdf("pdf/DLPFC_CellType_PTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "DLPFC PTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_dACC[,1]
mat <- data.matrix(li_dACC[,2:ncol(li_dACC)])
rownames(mat) <- rnames

pdf("pdf/dACC_CellType_PTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "dACC PTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device



########################################################################



li_BasoAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_onlyPTSD_Li2018.csv")
li_MedialAmyg <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_onlyPTSD_Li2018.csv")
li_dACC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_Li2018.csv")
li_DLPFC <-read.csv("csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_onlyPTSD_Li2018.csv")

colnames(li_DLPFC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_dACC) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_MedialAmyg) <-c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")
colnames(li_BasoAmyg) <- c("Type", "adj P 0.05", "adj P 0.01", "adj P 0.001", "adj P 0.0001")

#onlyPTSD
#Li
rnames <-li_BasoAmyg[,1]
mat <- data.matrix(li_BasoAmyg[,2:ncol(li_BasoAmyg)])
rownames(mat) <- rnames

pdf("pdf/BLA_CellType_onlyPTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "BLA onlyPTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_MedialAmyg[,1]
mat <- data.matrix(li_MedialAmyg[,2:ncol(li_MedialAmyg)])
rownames(mat) <- rnames

pdf("pdf/MedialAmyg_CellType_onlyPTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "MeA onlyPTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_DLPFC[,1]
mat <- data.matrix(li_DLPFC[,2:ncol(li_DLPFC)])
rownames(mat) <- rnames

pdf("pdf/DLPFC_CellType_onlyPTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "DLPFC onlyPTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device


rnames <-li_dACC[,1]
mat <- data.matrix(li_dACC[,2:ncol(li_dACC)])
rownames(mat) <- rnames

pdf("pdf/dACC_CellType_onlyPTSD_heatmap.pdf")
heatmap.2(mat,
  cellnote = round(mat,3),  # same data set for cell labels
  main = "dACC onlyPTSD", # heat map title
  notecol="black",      # change font color of cell labels to black
  key=FALSE,
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9))     # widens margins around plot
dev.off()               # close the PNG device

