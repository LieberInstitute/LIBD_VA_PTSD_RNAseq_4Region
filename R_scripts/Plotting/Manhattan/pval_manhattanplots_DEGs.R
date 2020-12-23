#Manhattan plots
#helpful other code to check out: https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function

#gist of what we're trying to get our object looking like:
#make df of DE analyses results 
##GENE CHR START END BP REGION +/-adj.P.Val (tab-delim)

#    unique gene identifier (gene_id, external_gene_name etc.) --> Symbol and ensemblID
#    "chr", "start", and "end"
#    chr length    
#    midpoint of gene
#    new start and end (additive for plotting)  
#    new midpoint of gene
#    brain region
#    the probability measure of significance (adj.P.Val)
#    logFC
#    final adj.P.Val (-log10(adj.P.Val),neg if neg logFC)


###################################################################
############PTSD
###################################################################

library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load rse object and get into form that was used for DE analysis
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#extract phenotype data
pd = colData(rse_gene)
#extract rowRanges
rr = data.frame(rowRanges(rse_gene))

#load DE results
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda",verbose=TRUE)
BLA_PTSD <- data.frame(geneStats_BasoAmygall[,grep("_PTSD\\b",colnames(geneStats_BasoAmygall))])
colnames(BLA_PTSD) <- paste0(colnames(BLA_PTSD),"_","BLA")

load("rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda",verbose=TRUE)
MeA_PTSD <- data.frame(geneStats_MedialAmygall[,grep("_PTSD\\b",colnames(geneStats_MedialAmygall))])
colnames(MeA_PTSD) <- paste0(colnames(MeA_PTSD),"_","MeA")

load("rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda",verbose=TRUE)
dACC_PTSD <- data.frame(geneStats_dACCall[,grep("_PTSD\\b",colnames(geneStats_dACCall))])
colnames(dACC_PTSD) <- paste0(colnames(dACC_PTSD),"_","dACC")

load("rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda",verbose=TRUE)
DLPFC_PTSD <- data.frame(geneStats_DLPFCall[,grep("_PTSD\\b",colnames(geneStats_DLPFCall))])
colnames(DLPFC_PTSD) <- paste0(colnames(DLPFC_PTSD),"_","DLPFC")

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda",verbose=TRUE)
All_PTSD <- data.frame(geneStatsall[,grep("_PTSD\\b",colnames(geneStatsall))])
colnames(All_PTSD) <- paste0(colnames(All_PTSD),"_","All")

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",verbose=TRUE)
Cortex_PTSD <- data.frame(geneStatsall[,grep("_PTSD\\b",colnames(geneStatsall))])
colnames(Cortex_PTSD) <- paste0(colnames(Cortex_PTSD),"_","Cortex")

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",verbose=TRUE)
Amyg_PTSD <- data.frame(geneStatsall[,grep("_PTSD\\b",colnames(geneStatsall))])
colnames(Amyg_PTSD) <- paste0(colnames(Amyg_PTSD),"_","Amyg")

seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:25]

#start combining
genes2match <- rownames(BLA_PTSD)
together<- cbind(BLA_PTSD, MeA_PTSD[match(genes2match,rownames(MeA_PTSD)),],dACC_PTSD[match(genes2match,rownames(dACC_PTSD)),],DLPFC_PTSD[match(genes2match,rownames(DLPFC_PTSD)),],
	All_PTSD[match(genes2match,rownames(All_PTSD)),], Cortex_PTSD[match(genes2match,rownames(Cortex_PTSD)),], Amyg_PTSD[match(genes2match,rownames(Amyg_PTSD)),])
testing<-cbind(BLA_PTSD, MeA_PTSD, dACC_PTSD,DLPFC_PTSD,All_PTSD,Cortex_PTSD,Amyg_PTSD)
identical(together, testing, attrib.as.set=FALSE)
#true

identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_MeA,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_dACC,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_DLPFC,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_All,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_Cortex,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_Amyg,attrib.as.set=FALSE)
#[1] TRUE

together_less <- together

together_less$Symbol <- rr$Symbol[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$ensemblID <- rr$ensemblID[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$chr <- rr$seqnames[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$start <- rr$start[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$end <- rr$end[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$width <- rr$width[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$gene_length <- rr$Length[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
building <- together_less

building$running <- (building$start + building$end)/2

chr_vals <-seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:25]
#Note, there are some unlocalized segments in the TxDb, but I dont think we need to incorporate that info for these plots
chr_vals<- as.numeric(chr_vals)

chr1 <- 0
chr2 <- chr_vals[1]
chr3 <- chr_vals[1] + chr_vals[2]
chr4 <- chr_vals[1] + chr_vals[2] + chr_vals[3]
chr5 <- chr4 + chr_vals[4]
chr6 <- chr5 + chr_vals[5]
chr7 <- chr6 + chr_vals[6]
chr8 <- chr7 + chr_vals[7]
chr9 <- chr8 + chr_vals[8]
chr10 <- chr9 + chr_vals[9]
chr11 <- chr10 + chr_vals[10]
chr12 <- chr11 + chr_vals[11]
chr13 <- chr12 + chr_vals[12]
chr14 <- chr13 + chr_vals[13]
chr15 <- chr14 + chr_vals[14]
chr16 <- chr15 + chr_vals[15]
chr17 <- chr16 + chr_vals[16]
chr18 <- chr17 + chr_vals[17]
chr19 <- chr18 + chr_vals[18]
chr20 <- chr19 + chr_vals[19]
chr21 <- chr20 + chr_vals[20]
chr22 <- chr21 + chr_vals[21]
chrX <- chr22 + chr_vals[22]
chrY <- chrX + chr_vals[23]
chrM <- chrY + chr_vals[24]

all_chr_vals <- c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM)

names(all_chr_vals) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
	"chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
	"chrX","chrY","chrM")

building$chr <- as.character(building$chr)
building$all_chr_vals <- sapply(building$chr, function(x) all_chr_vals[which(names(all_chr_vals) %in% x)])
building$gene_loc <- building$all_chr_vals + building$running

building$final.P.adj_BLA <- (-log10(building$adj.P.Val_PTSD_BLA))
building$final.P.adj_MeA <- (-log10(building$adj.P.Val_PTSD_MeA))
building$final.P.adj_dACC <- (-log10(building$adj.P.Val_PTSD_dACC))
building$final.P.adj_DLPFC <- (-log10(building$adj.P.Val_PTSD_DLPFC))
building$final.P.adj_All <- (-log10(building$adj.P.Val_PTSD_All))
building$final.P.adj_Cortex <- (-log10(building$adj.P.Val_PTSD_Cortex))
building$final.P.adj_Amyg <- (-log10(building$adj.P.Val_PTSD_Amyg))

building$final.P.adj_BLA <- ifelse(building$logFC_PTSD_BLA > 0, building$final.P.adj_BLA, -building$final.P.adj_BLA)
building$final.P.adj_MeA <- ifelse(building$logFC_PTSD_MeA > 0, building$final.P.adj_MeA, -building$final.P.adj_MeA)
building$final.P.adj_dACC <- ifelse(building$logFC_PTSD_dACC > 0, building$final.P.adj_dACC, -building$final.P.adj_dACC)
building$final.P.adj_DLPFC <- ifelse(building$logFC_PTSD_DLPFC > 0, building$final.P.adj_DLPFC, -building$final.P.adj_DLPFC)
building$final.P.adj_All <- ifelse(building$logFC_PTSD_All > 0, building$final.P.adj_All, -building$final.P.adj_All)
building$final.P.adj_Cortex <- ifelse(building$logFC_PTSD_Cortex > 0, building$final.P.adj_Cortex, -building$final.P.adj_Cortex)
building$final.P.adj_Amyg <- ifelse(building$logFC_PTSD_Amyg > 0, building$final.P.adj_Amyg, -building$final.P.adj_Amyg)

building$final.P.Value_BLA <- (-log10(building$P.Value_PTSD_BLA))
building$final.P.Value_MeA <- (-log10(building$P.Value_PTSD_MeA))
building$final.P.Value_dACC <- (-log10(building$P.Value_PTSD_dACC))
building$final.P.Value_DLPFC <- (-log10(building$P.Value_PTSD_DLPFC))
building$final.P.Value_All <- (-log10(building$P.Value_PTSD_All))
building$final.P.Value_Cortex <- (-log10(building$P.Value_PTSD_Cortex))
building$final.P.Value_Amyg <- (-log10(building$P.Value_PTSD_Amyg))

building$final.P.Value_BLA <- ifelse(building$logFC_PTSD_BLA > 0, building$final.P.Value_BLA, -building$final.P.Value_BLA)
building$final.P.Value_MeA <- ifelse(building$logFC_PTSD_MeA > 0, building$final.P.Value_MeA, -building$final.P.Value_MeA)
building$final.P.Value_dACC <- ifelse(building$logFC_PTSD_dACC > 0, building$final.P.Value_dACC, -building$final.P.Value_dACC)
building$final.P.Value_DLPFC <- ifelse(building$logFC_PTSD_DLPFC > 0, building$final.P.Value_DLPFC, -building$final.P.Value_DLPFC)
building$final.P.Value_All <- ifelse(building$logFC_PTSD_All > 0, building$final.P.Value_All, -building$final.P.Value_All)
building$final.P.Value_Cortex <- ifelse(building$logFC_PTSD_Cortex > 0, building$final.P.Value_Cortex, -building$final.P.Value_Cortex)
building$final.P.Value_Amyg <- ifelse(building$logFC_PTSD_Amyg > 0, building$final.P.Value_Amyg, -building$final.P.Value_Amyg)


###########################################################################################################


length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_BLA < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_MeA < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_dACC < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_DLPFC < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_All < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_Cortex < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_Amyg < 0.3])
#0
dim(building)
#26020   136

building <- building[!building$chr == "chrM",]
dim(building)
#25983   136


##color based on adjusted P
#FDR < 0.1, special color
#FDR < 0.05, special shape (if TRUE, triangle. if FALSE, circle)
#if pos logFC and triangle: triangle points up, if neg logFC and triangle: triangle points down

#BLA
building$BLA_color <- ifelse(building$final.P.adj_BLA > 1, "#228B22",ifelse(building$final.P.adj_BLA < -1, "#228B22",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#9ECAE1",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#2171B5",NA))))

#MeA
building$MeA_color <- ifelse(building$final.P.adj_MeA > 1, "#EE7600",ifelse(building$final.P.adj_MeA < -1, "#EE7600",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

#dACC
building$dACC_color <- ifelse(building$final.P.adj_dACC > 1, "#FF34B3",ifelse(building$final.P.adj_dACC < -1, "#FF34B3",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

#DLPFC
building$DLPFC_color <- ifelse(building$final.P.adj_DLPFC > 1, "#B23AEE",ifelse(building$final.P.adj_DLPFC < -1, "#B23AEE",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

#All
building$All_color <- ifelse(building$final.P.adj_All > 1, "#CD2626",ifelse(building$final.P.adj_All < -1, "#CD2626",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

#Cortex
building$Cortex_color <- ifelse(building$final.P.adj_Cortex > 1, "#C43EB4",ifelse(building$final.P.adj_Cortex < -1, "#C43EB4",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

#Amyg
building$Amyg_color <- ifelse(building$final.P.adj_Amyg > 1, "#FFCF00",ifelse(building$final.P.adj_Amyg < -1, "#FFCF00",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

bla_col_values <-c("#228B22" = "#228B22", "#2171B5" = "#2171B5","#9ECAE1" = "#9ECAE1")
mea_col_values <-c("#EE7600" = "#EE7600", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
dacc_col_values <-c("#FF34B3" = "#FF34B3", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
dlpfc_col_values <-c("#B23AEE" = "#B23AEE", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
All_col_values <-c("#CD2626" = "#CD2626", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
Cortex_col_values <-c("#C43EB4" = "#C43EB4", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
Amyg_col_values <-c("#FFCF00" = "#FFCF00", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")

#region_palette:
#CD2626: All (across, PTSD term)
#228B22: BLA
#EE7600: MeA
#B23AEE: DLPFC
#FF34B3: dACC
#F0C900 : amyg
#C43EB4 : cortex

#triangle: pch = 17,24,25
#circle: pch = 19,21

font="Helvetica"
text_color="#222222"

#get x axis label positions
chr_vals <-seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:25]
chr_vals<- as.numeric(chr_vals)
x_vals <- chr_vals[1:24]

x1 <- 0.5*x_vals[1]
x2 <- x_vals[1] + 0.5*x_vals[2]
x3 <- x_vals[1] + x_vals[2] + 0.5*x_vals[3]
x4 <- x_vals[1] + x_vals[2] + x_vals[3] + 0.5*x_vals[4]
x5 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + 0.5*x_vals[5]
x6 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + 0.5*x_vals[6]
x7 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + 0.5*x_vals[7]
x8 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + 0.5*x_vals[8]
x9 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + 0.5*x_vals[9]
x10 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + 0.5*x_vals[10]
x11 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + 0.5*x_vals[11]
x12 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] +0.5*x_vals[12]
x13 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12] +0.5*x_vals[13]
x14 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13] +0.5*x_vals[14]
x15 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14] +0.5*x_vals[15]
x16 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15] +0.5*x_vals[16]
x17 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16] +0.5*x_vals[17]
x18 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + 0.5*x_vals[18]
x19 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] +0.5*x_vals[19]
x20 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] + x_vals[19] +0.5*x_vals[20]
x21 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] + x_vals[19] + x_vals[20] +0.5*x_vals[21]
x22 <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] + x_vals[19] + x_vals[20] + x_vals[21] +0.5*x_vals[22]
xX <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] + x_vals[19] + x_vals[20] + x_vals[21] + x_vals[22] +0.5*x_vals[23]
xY <- x_vals[1] + x_vals[2] + x_vals[3] + x_vals[4] + x_vals[5] + x_vals[6] + x_vals[7] + x_vals[8] + x_vals[9] + x_vals[10] + x_vals[11] + x_vals[12]+ x_vals[13]+ x_vals[14]+ x_vals[15]+ x_vals[16]+ x_vals[17] + x_vals[18] + x_vals[19] + x_vals[20] + x_vals[21] + x_vals[22] + x_vals[23] +0.5*max(building$end[building$chr == "chrY"])
#note, for xY, i used half the max of our farthest reaching gene bc there seemed to be no gene info for a lot of the remaining length of chrY based on x_vals[24]
all_x_vals <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,xX,xY)


#FDR < 0.1 symbols found in amyg, cortex, and allregions

intersect(building$Symbol_PTSD_All[building$adj.P.Val_PTSD_All < 0.1],
	intersect(building$Symbol_PTSD_Cortex[building$adj.P.Val_PTSD_Cortex < 0.1],building$Symbol_PTSD_Amyg[building$adj.P.Val_PTSD_Amyg < 0.1]))
#none

intersect(building$Symbol_PTSD_Cortex[building$adj.P.Val_PTSD_Cortex < 0.1],building$Symbol_PTSD_Amyg[building$adj.P.Val_PTSD_Amyg < 0.1])
#none

all_amyg <- intersect(building$Symbol_PTSD_All[building$adj.P.Val_PTSD_All < 0.1],building$Symbol_PTSD_Amyg[building$adj.P.Val_PTSD_Amyg < 0.1])
#[1] "CRHBP"

all_cortex <- intersect(building$Symbol_PTSD_All[building$adj.P.Val_PTSD_All < 0.1],building$Symbol_PTSD_Cortex[building$adj.P.Val_PTSD_Cortex < 0.1])
#33 genes

allsyms <- c(all_amyg, all_cortex)

save(building, all_chr_vals, all_x_vals, allsyms,bla_col_values,mea_col_values,dacc_col_values,dlpfc_col_values,All_col_values,Cortex_col_values,Amyg_col_values,
 font, text_color, file="rdas/Manhattan_plots_Pvals_PTSD_objects.rda")


load("rdas/Manhattan_plots_Pvals_PTSD_objects.rda",verbose=TRUE)

#####################################
#BasoAmyg
#####################################
##color and shape based on adjusted P
#FDR < 0.1, special color
#FDR < 0.05, special shape (if TRUE, triangle. if FALSE, circle)
#if pos logFC and triangle: up, if neg logFC and triangle: down

##hlines based on adjusted P such that i found the p value corresponding to FDR = 0.1 (approximately), calc for each region separately
#P = 1.568715e-05 for FDR = 0.1 (-log10(P) = 4.804456)
#calclated the above by taking the values immediately around adjusted P = 0.1 (since 0.1 wasnt exactly available from the data). I averaged the corresponding p values if the 
#immediate FDR values were non-unique.Then I used y2-y0 = slope(x2-x0) to calc y2 (new p value) given known x (0.1). slope was calc as (y1-y0)/(x1-x0)

manhattan_BLA <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_BLA, x=building$gene_loc, color=building$BLA_color), 
shape=ifelse(building$adj.P.Val_PTSD_BLA < 0.05 & building$final.P.adj_BLA < 0, 25, ifelse(building$adj.P.Val_PTSD_BLA < 0.05 & building$final.P.adj_BLA > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_BLA < 0.05, 6,4), fill=building$BLA_color) +
scale_color_manual(values=bla_col_values) +
scale_fill_manual(values=bla_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=4.804456, color="#228B22",size=0.75) +
geom_hline(yintercept=-4.804456,color="#228B22",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = subset(building, adj.P.Val_PTSD_BLA <0.1), aes(label = Symbol, y=final.P.Value_BLA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_x=50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Basolateral Amygdala") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_BLA,file="png/Manhattan_BasoAmyg_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


#####################################
#MedialAmyg
#####################################
#FDR = 0.1 --> P = 6.17859e-05, -log10P = 4.209111

manhattan_MeA <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_MeA, x=building$gene_loc, color=building$MeA_color),
shape=ifelse(building$adj.P.Val_PTSD_MeA < 0.05 & building$final.P.adj_MeA < 0, 25, ifelse(building$adj.P.Val_PTSD_MeA < 0.05 & building$final.P.adj_MeA > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_MeA < 0.05, 6,4), fill=building$MeA_color) +
scale_color_manual(values=mea_col_values) +
scale_fill_manual(values=mea_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=4.209111, color="#EE7600",size=0.75) +
geom_hline(yintercept=-4.209111,color="#EE7600",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_MeA,decreasing=TRUE),][c(1:5,7:8),], aes(label = Symbol, y=final.P.Value_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_MeA,decreasing=TRUE),][6,], aes(label = Symbol, y=final.P.Value_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),vjust=1,hjust=0,nudge_y=0.75) +
geom_text_repel(data = subset(building, final.P.adj_MeA < -1)[c(1:5,8),], aes(label = Symbol, y=final.P.Value_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_y=-0.75) +
geom_text_repel(data = subset(building, final.P.adj_MeA < -1)[c(6:7),], aes(label = Symbol, y=final.P.Value_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_x = 50000000,hjust=0,vjust=0.75) +
ggtitle("Differentially Expressed Genes", subtitle="Medial Amygdala") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_MeA,file="png/Manhattan_MedialAmyg_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


#####################################
#dACC
#####################################
##hlines based on adjusted P such that i found the p value corresponding to FDR = 0.1 (approximately), calc for each region separately
#P = 0.0002751387 for FDR = 0.1 (-log10(P) = 3.560448)
#calclated the above by taking the values immediately around adjusted P = 0.1 (since 0.1 wasnt exactly available from the data). I averaged the corresponding p values if the 
#immediate FDR values were non-unique.Then I used y2-y0 = slope(x2-x0) to calc y2 (new p value) given known x (0.1). slope was calc as (y1-y0)/(x1-x0)

manhattan_dACC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_dACC, x=building$gene_loc, color=building$dACC_color),
shape=ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC < 0, 25, ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_dACC < 0.05, 6,4), fill=building$dACC_color) +
scale_color_manual(values=dacc_col_values) +
scale_fill_manual(values=dacc_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=3.560448, color="#FF34B3",size=0.75) +
geom_hline(yintercept=-3.560448,color="#FF34B3",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(1),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(7),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = 50000000,hjust=0,vjust=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(5),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.70,vjust=1,hjust=0) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(2,3,4,9,10),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(6,8),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.80) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(2,5,7,8,10,11),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(1,3,4,6,9),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = 50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Dorsal Anterior Cingulate Cortex") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_dACC,file="png/Manhattan_dACC_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


##also include FDR < 0.05 line
#P is approx 4.048699e-05, -log10 = 4.392685

manhattan_dACC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_dACC, x=building$gene_loc, color=building$dACC_color),
shape=ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC < 0, 25, ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC > 0, 24,21)),size=4, 
fill=building$dACC_color) +
geom_point(aes(y=rep(4.392685,103), x=seq(0,3080000000,by=30000000)), color="#222222",shape=24,size=1.75, fill="#222222",alpha=0.75) +
geom_point(aes(y=rep(-4.392685,103), x=seq(0,3080000000,by=30000000)), color="#222222",shape=25,size=1.75, fill="#222222",alpha=0.75) +
scale_color_manual(values=dacc_col_values) +
scale_fill_manual(values=dacc_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=3.560448, color="#FF34B3",size=.75) +
geom_hline(yintercept=-3.560448,color="#FF34B3",size=.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(1),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(7),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.75,hjust=0,vjust=1) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(5),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.75,vjust=1,hjust=0) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(2,3,4,9,10),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(8),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(6),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.8,hjust=1,vjust=1) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(2,5,7,8,10),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.75) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(1,3,4,6),], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][9,], aes(label = Symbol, y=final.P.Value_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-1.1) +
ggtitle("Differentially Expressed Genes", subtitle="Dorsal Anterior Cingulate Cortex") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_dACC,file="png/Manhattan_dACC_PTSD_labels_Pvals_twohlines.png",width=16,height=12,dpi=320)

########################################################################
#compare to other plot (previous) type
manhattan_dACC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_dACC, x=building$gene_loc, color=building$dACC_color),
shape=ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC < 0, 25, ifelse(building$adj.P.Val_PTSD_dACC < 0.05 & building$final.P.adj_dACC > 0, 24,21)),size=4, 
fill=building$dACC_color) +
scale_color_manual(values=dacc_col_values) +
scale_fill_manual(values=dacc_col_values) +
scale_y_continuous(limits=c(-2.75,2.75),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-2.5, 2.5, by = 0.5)) +
geom_hline(yintercept=1, color="#FF34B3",size=0.75) +
geom_hline(yintercept=-1,color="#FF34B3",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(1),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.25) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(7),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(5),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.20,vjust=1,hjust=0) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(2,3,4,9,10),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.25) +
geom_text_repel(data = building[order(building$final.P.adj_dACC,decreasing=TRUE),][c(6,8),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.3) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(2,5,7,8,10),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.25) +
geom_text_repel(data = building[order(building$final.P.adj_dACC),][c(1,3,4,6,9),], aes(label = Symbol, y=final.P.adj_dACC, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = 50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Dorsal Anterior Cingulate Cortex") + 
ylab("-log10(Adjusted P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_dACC,file="png/Manhattan_dACC_PTSD_labels_adjP.png",width=16,height=12,dpi=320)
################################################################################################################


#####################################
#DLPFC
#####################################
#FDR = 0.1, p = 1.379053e-05, -log10(p) = 4.860419

manhattan_DLPFC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_DLPFC, x=building$gene_loc, color=building$DLPFC_color),
shape=ifelse(building$adj.P.Val_PTSD_DLPFC < 0.05 & building$final.P.adj_DLPFC < 0, 25, ifelse(building$adj.P.Val_PTSD_DLPFC < 0.05 & building$final.P.adj_DLPFC > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_DLPFC < 0.05, 6,4), fill=building$DLPFC_color) +
scale_color_manual(values=dlpfc_col_values) +
scale_fill_manual(values=dlpfc_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=4.860419, color="#B23AEE",size=0.75) +
geom_hline(yintercept=-4.860419,color="#B23AEE",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = subset(building, adj.P.Val_PTSD_DLPFC <0.1), aes(label = Symbol, y=final.P.Value_DLPFC, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_x=50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Dorsolateral Prefrontal Cortex") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_DLPFC,file="png/Manhattan_DLPFC_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


#####################################
#All Regions
#####################################
#FDR = 0.1, p = 0.001061907, -log10P = 2.973914
keepsyms <-allsyms[c(1,2,7,22)]

intersect(keepsyms,building$Symbol[order(building$final.P.adj_All,decreasing=TRUE)][c(1:10)])
#none
intersect(keepsyms,building$Symbol[order(building$final.P.adj_All)][c(1:10)])
#CRHBP

diffsyms <- setdiff(keepsyms,building$Symbol[order(building$final.P.adj_All,decreasing=TRUE)][c(1:10)])
diffsyms <- setdiff(diffsyms,building$Symbol[order(building$final.P.adj_All)][c(1:10)])

manhattan_All <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_All, x=building$gene_loc, color=building$All_color),
shape=ifelse(building$adj.P.Val_PTSD_All < 0.05 & building$final.P.adj_All < 0, 25, ifelse(building$adj.P.Val_PTSD_All < 0.05 & building$final.P.adj_All > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_All < 0.05, 6,4), fill=building$All_color) +
scale_color_manual(values=All_col_values) +
scale_fill_manual(values=All_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=2.973914, color="#CD2626",size=1) +
geom_hline(yintercept=-2.973914,color="#CD2626",size=1) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +

geom_text_repel(data = building[order(building$final.P.adj_All,decreasing=TRUE),][c(1,3,4,8,9),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_All,decreasing=TRUE),][c(7),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_All,decreasing=TRUE),][c(6),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_All,decreasing=TRUE),][c(2,5,10),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=-0.75) +
#geom_text_repel(data = building[building$Symbol %in% diffsyms & building$final.P.Value_All > 0,], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),
#	min.segment.length = unit(0, 'lines'),color=text_color, family=font,fontface='bold',size=8,nudge_y=0.75) +
geom_text_repel(data = building[order(building$final.P.adj_All),][c(2:8,10),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=-0.75) +
geom_text_repel(data = building[order(building$final.P.adj_All),][c(1,9),], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$Symbol %in% "FOS" & building$final.P.Value_All < 0,], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),
	min.segment.length = unit(0, 'lines'),color=text_color, family=font,fontface='bold',size=8,nudge_y=-0.65) +
geom_text_repel(data = building[building$Symbol %in% c("CORT") & building$final.P.Value_All < 0,], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),
	min.segment.length = unit(0, 'lines'),color=text_color, family=font,fontface='bold',size=8,nudge_x = 150000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$Symbol %in% c("SST") & building$final.P.Value_All < 0,], aes(label = Symbol, y=final.P.Value_All, x=gene_loc),
	min.segment.length = unit(0, 'lines'),color=text_color, family=font,fontface='bold',size=8,nudge_x = -70000000,hjust=1,vjust=0.5) +

ggtitle("Differentially Expressed Genes", subtitle="Four Region Interaction Model, PTSD Term") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_All,file="png/Manhattan_4RegionIntxn_GroupPTSD_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


#####################################
#Cortex
#####################################
#FDR = 0.1, p =0.0004576275, -log10P = 3.339488

keepsyms <-allsyms[c(1,2,7,22)]

intersect(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.adj_Cortex < 0])
#SST and CORT
intersect(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.adj_Cortex > 0])
#none

diffsyms <- setdiff(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Cortex < 0.05])

manhattan_Cortex <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_Cortex, x=building$gene_loc, color=building$Cortex_color),
shape=ifelse(building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.adj_Cortex < 0, 25, ifelse(building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.adj_Cortex > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_Cortex < 0.05, 6,4), fill=building$Cortex_color) +
scale_color_manual(values=Cortex_col_values) +
scale_fill_manual(values=Cortex_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=3.339488, color="#C43EB4",size=1) +
geom_hline(yintercept=-3.339488,color="#C43EB4",size=1) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +

geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(1,4:8,16,18),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=0.75) + #18
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(3),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(2),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=2) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(14),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=1.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(9),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = -50000000,hjust=1,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(11),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = -100000000,hjust=1,vjust=0) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(13),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=2.4) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(12),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=1.2,hjust=1,vjust=0,nudge_x = -100000000) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(15),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=2.75) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(17),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex > 0,][c(10),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=0.5,hjust=1,vjust=0,nudge_x = -100000000) +

geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(1,2),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -0.75) + #23
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(3),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.2,hjust=0,vjust=1,nudge_x=200000000) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(4),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.25) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(5),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=-0.50,hjust=0,vjust=0.5,nudge_x=50000000) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(6),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = -40000000,hjust=1,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(8),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -2.25) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(9),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -2.75) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(10),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1,hjust=0,vjust=1,nudge_x=100000000) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(13),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.15) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(14),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -2) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(18),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = -70000000,hjust=1,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(22,23),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(21),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.75,hjust=0,vjust=1,nudge_x = 220000000) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(19),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -2.25) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(20),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.25) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(12),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 50000000,hjust=0,vjust=0.5) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(17),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -1.3) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Cortex < 0.05 & building$final.P.Value_Cortex < 0,][c(7,11,15,16),], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -0.75) +

#geom_text_repel(data = building[building$Symbol %in% diffsyms & building$final.P.Value_Cortex > 0,], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
#	family=font,fontface='bold',size=8,nudge_y = 0.75) +

geom_text_repel(data = building[building$Symbol %in% "CRHBP" & building$final.P.Value_Cortex < 0,], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -0.75,hjust=1,vjust=1,nudge_x=-70000000) +
geom_text_repel(data = building[building$Symbol %in% "FOS" & building$final.P.Value_Cortex < 0,], aes(label = Symbol, y=final.P.Value_Cortex, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = -50000000,hjust=1,vjust=1) +

ggtitle("Differentially Expressed Genes", subtitle="Four Region Interaction Model, Frontal Cortex, PTSD Term") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_Cortex,file="png/Manhattan_4RegionIntxn_Cortex_GroupPTSD_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)


#####################################
#Amyg
#####################################
#FDR = 0.1, p =7.669214e-06, -log10P = 5.115249

keepsyms <-allsyms[c(1,2,7,22)]

intersect(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.adj_Amyg < 0])
#CRHBP
intersect(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.adj_Amyg > 0])
#none

diffsyms <- setdiff(keepsyms,building$Symbol[building$adj.P.Val_PTSD_Amyg < 0.05])


manhattan_Amyg <- ggplot() +
geom_point(data=building,aes(y=building$final.P.Value_Amyg, x=building$gene_loc, color=building$Amyg_color),
shape=ifelse(building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.adj_Amyg < 0, 25, ifelse(building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.adj_Amyg > 0, 24,21)),
size=ifelse(building$adj.P.Val_PTSD_Amyg < 0.05, 6,4), fill=building$Amyg_color) +
scale_color_manual(values=Amyg_col_values) +
scale_fill_manual(values=Amyg_col_values) +
scale_y_continuous(limits=c(-8,8),labels = c(seq(8, 0, by = -1),seq(1, 8, by = 1)), breaks=seq(-8, 8, by = 1)) +
geom_hline(yintercept=5.115249, color="#FFCF00",size=1) +
geom_hline(yintercept=-5.115249,color="#FFCF00",size=1) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +

geom_text_repel(data = building[building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.Value_Amyg > 0,], aes(label = Symbol, y=final.P.Value_Amyg, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y=0.75) +
geom_text_repel(data = building[building$adj.P.Val_PTSD_Amyg < 0.05 & building$final.P.Value_Amyg < 0,], aes(label = Symbol, y=final.P.Value_Amyg, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_x = 70000000,hjust=0,vjust=0.5) +
#geom_text_repel(data = building[building$Symbol %in% diffsyms & building$final.P.Value_Amyg > 0,], aes(label = Symbol, y=final.P.Value_Amyg, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
#	family=font,fontface='bold',size=8,nudge_y = 0.75) +
geom_text_repel(data = building[building$Symbol %in% diffsyms & building$final.P.Value_Amyg < 0,], aes(label = Symbol, y=final.P.Value_Amyg, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=8,nudge_y = -0.75) +

ggtitle("Differentially Expressed Genes", subtitle="Four Region Interaction Model, Amygdala, PTSD Term") + 
ylab("-log10(P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = c(seq(1, 22, by = 1),"X","Y"), breaks=all_x_vals) +
theme(
	panel.background = element_blank(),
	panel.grid.major.y = element_line(colour="lightgrey",linetype="dashed"),
	panel.grid.major.x = element_blank(),
	panel.grid.minor = element_blank(), 
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(t=0.25,r=0,b=5,l=0),hjust = 0.5),
	legend.position = "none",
	#legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(t=15)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin=ggplot2::margin(r=15)),
	axis.text.x = element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(t=5)),
	axis.ticks.x=element_line(colour = text_color,size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_text(family=font, size=16,color=text_color, face="bold",margin=ggplot2::margin(r=5)),
	axis.ticks.y=element_blank())
ggsave(plot=manhattan_Amyg,file="png/Manhattan_4RegionIntxn_Amyg_GroupPTSD_PTSD_labels_Pvals.png",width=16,height=12,dpi=320)

