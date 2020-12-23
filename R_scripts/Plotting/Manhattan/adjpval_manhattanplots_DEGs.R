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

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_genealone_threeGroup.rda",verbose=TRUE)
Pan_PTSD <- data.frame(geneStatsall[,grep("_PTSD\\b",colnames(geneStatsall))])
colnames(Pan_PTSD) <- paste0(colnames(Pan_PTSD),"_","Pan")

seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[1:25]

#start combining
genes2match <- rownames(BLA_PTSD)
together<- cbind(BLA_PTSD, MeA_PTSD[match(genes2match,rownames(MeA_PTSD)),],dACC_PTSD[match(genes2match,rownames(dACC_PTSD)),],DLPFC_PTSD[match(genes2match,rownames(DLPFC_PTSD)),],
	Pan_PTSD[match(genes2match,rownames(Pan_PTSD)),])
testing<-cbind(BLA_PTSD, MeA_PTSD, dACC_PTSD,DLPFC_PTSD,Pan_PTSD)
identical(together, testing, attrib.as.set=FALSE)
#true

identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_MeA,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_dACC,attrib.as.set=FALSE)
#[1] TRUE
identical(together$ensemblID_PTSD_BLA,together$ensemblID_PTSD_DLPFC,attrib.as.set=FALSE)
#[1] TRUE

together_subset <- together[,c(1:7,11,15,27,31,43,47,59,63,75,79)]
together_less <- together_subset[,c(1,3,8:17)]

together_less$Symbol <- rr$Symbol[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$ensemblID <- rr$ensemblID[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$chr <- rr$seqnames[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$start <- rr$start[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$end <- rr$end[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$width <- rr$width[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
together_less$gene_length <- rr$Length[match(together_less$ensemblID_PTSD_BLA, rr$ensemblID)]
building <- together_less[,c(13:19,3:12)]

colnames(building)
# [1] "Symbol"               "ensemblID"            "chr"                 
# [4] "start"                "end"                  "width"               
# [7] "gene_length"          "logFC_PTSD_BLA"       "adj.P.Val_PTSD_BLA"  
#[10] "logFC_PTSD_MeA"       "adj.P.Val_PTSD_MeA"   "logFC_PTSD_dACC"     
#[13] "adj.P.Val_PTSD_dACC"  "logFC_PTSD_DLPFC"     "adj.P.Val_PTSD_DLPFC"
#[16] "logFC_PTSD_Pan"       "adj.P.Val_PTSD_Pan" 

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
building$final.P.adj_Pan <- (-log10(building$adj.P.Val_PTSD_Pan))

building$final.P.adj_BLA <- ifelse(building$logFC_PTSD_BLA > 0, building$final.P.adj_BLA, -building$final.P.adj_BLA)
building$final.P.adj_MeA <- ifelse(building$logFC_PTSD_MeA > 0, building$final.P.adj_MeA, -building$final.P.adj_MeA)
building$final.P.adj_dACC <- ifelse(building$logFC_PTSD_dACC > 0, building$final.P.adj_dACC, -building$final.P.adj_dACC)
building$final.P.adj_DLPFC <- ifelse(building$logFC_PTSD_DLPFC > 0, building$final.P.adj_DLPFC, -building$final.P.adj_DLPFC)
building$final.P.adj_Pan <- ifelse(building$logFC_PTSD_Pan > 0, building$final.P.adj_Pan, -building$final.P.adj_Pan)
								   
save(building, all_chr_vals, file="rdas/Manhattan_plots_adjPvals_PTSD_objects.rda")

###########################################################################################################

load("rdas/Manhattan_plots_adjPvals_PTSD_objects.rda", verbose=TRUE)

length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_BLA < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_MeA < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_dACC < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_DLPFC < 0.3])
#0
length(building$Symbol[building$chr == "chrM" & building$adj.P.Val_PTSD_Pan < 0.3])
#1 --> MT-TT, adj P = 0.1443586

dim(building)
#[1] 26020    25

building <- building[!building$chr == "chrM",]
dim(building)
#[1] 25983    25

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

#Pan
building$Pan_color <- ifelse(building$final.P.adj_Pan > 1, "#CD2626",ifelse(building$final.P.adj_Pan < -1, "#CD2626",ifelse(building$chr %in%
	c("chr1","chr5","chr9","chr13","chr17","chr21","chr3","chr7","chr11","chr15","chr19","chrX"), "#2171B5",ifelse(building$chr %in%
	c("chr2","chr6","chr10","chr14","chr18","chr22","chr4","chr8","chr12","chr16","chr20","chrY"), "#9ECAE1",NA))))

bla_col_values <-c("#228B22" = "#228B22", "#2171B5" = "#2171B5","#9ECAE1" = "#9ECAE1")
mea_col_values <-c("#EE7600" = "#EE7600", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
dacc_col_values <-c("#FF34B3" = "#FF34B3", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
dlpfc_col_values <-c("#B23AEE" = "#B23AEE", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")
pan_col_values <-c("#CD2626" = "#CD2626", "#2171B5" = "#2171B5", "#9ECAE1" = "#9ECAE1")

#region_palette:
#CD2626: Pan (across, PTSD term)
#228B22: BLA
#EE7600: MeA
#B23AEE: DLPFC
##FF34B3: dACC

#triangle: pch = 17
#circle: pch = 19

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

#####################################
#BasoAmyg
#####################################

#with minor axis
hmm <- c(all_chr_vals, all_x_vals)
hmm<-as.numeric(hmm)
hmm<-sort(hmm)
labelz <- c("","1","", "2","", "3","", "4","", "5","", "6","", "7","", "8","" ,"9","" ,"10", "","11","", "12","", 
	"13","", "14","", "15","","16","", "17" ,"","18","","19" ,"","20","", "21" ,"","22","","X","","Y","")

manhattan_BLA_minor <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_BLA, x=building$gene_loc, color=building$BLA_color),size=3) +
scale_colour_manual(values=bla_col_values) +
scale_y_continuous(limits=c(-1.25,1.25),labels = c(seq(1.25, 0, by = -0.25),seq(0.25, 1.25, by = 0.25)), breaks=seq(-1.25, 1.25, by = 0.25)) +
geom_hline(yintercept=1, color="#228B22",size=0.75) +
geom_hline(yintercept=-1,color="#228B22",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = subset(building, adj.P.Val_PTSD_BLA <0.1), aes(label = Symbol, y=final.P.adj_BLA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6) +
ggtitle("Differentially Expressed Genes", subtitle="Basolateral Amygdala") + 
ylab("-log10(Adjusted P Value)") +
xlab("Position") +
scale_x_continuous(expand = c(0.001,0.001),limits=c(-30000000,3080000000),labels = labelz, breaks=hmm) +
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
ggsave(plot=manhattan_BLA_minor,file="png/Manhattan_BasoAmyg_PTSD_minoraxis_adjP.png",width=16,height=12,dpi=300)

#without minor axis
manhattan_BLA <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_BLA, x=building$gene_loc, color=building$BLA_color), shape=ifelse(building$final.P.adj_BLA < -1, 25, 
	ifelse(building$final.P.adj_BLA > 1, 24, 21)),size=ifelse(building$adj.P.Val_PTSD_BLA < 0.1, 4.5, 3.5), fill=building$BLA_color) +
scale_color_manual(values=bla_col_values) +
scale_fill_manual(values=bla_col_values) +
scale_y_continuous(limits=c(-1.25,1.25),labels = c(seq(1.25, 0, by = -0.25),seq(0.25, 1.25, by = 0.25)), breaks=seq(-1.25, 1.25, by = 0.25)) +
geom_hline(yintercept=1, color="#228B22",size=0.75) +
geom_hline(yintercept=-1,color="#228B22",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = subset(building, adj.P.Val_PTSD_BLA <0.1), aes(label = Symbol, y=final.P.adj_BLA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_x=50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Basolateral Amygdala") + 
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
ggsave(plot=manhattan_BLA,file="png/Manhattan_BasoAmyg_PTSD_labels_adjP.png",width=16,height=12,dpi=320)


#####################################
#MedialAmyg
#####################################

#without minor axis
manhattan_MeA <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_MeA, x=building$gene_loc, color=building$MeA_color),size=3) +
scale_colour_manual(values=mea_col_values) +
scale_y_continuous(limits=c(-1.25,1.25),labels = c(seq(1.25, 0, by = -0.25),seq(0.25, 1.25, by = 0.25)), breaks=seq(-1.25, 1.25, by = 0.25)) +
geom_hline(yintercept=1, color="#EE7600",size=0.75) +
geom_hline(yintercept=-1,color="#EE7600",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_MeA,decreasing=TRUE),][c(1:5,7:8),], aes(label = Symbol, y=final.P.adj_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_y=c(0.15,0.15,0.15,0.15,0.15,-0.15,0.15)) +
geom_text_repel(data = building[order(building$final.P.adj_MeA,decreasing=TRUE),][6,], aes(label = Symbol, y=final.P.adj_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),vjust=1,hjust=0,nudge_y=-0.1) +
geom_text_repel(data = subset(building, final.P.adj_MeA <= -1), aes(label = Symbol, y=final.P.adj_MeA, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_y=c(-0.15,0.15,-0.15,-0.15,0.15,-0.15,0.15,-0.15)) +
ggtitle("Differentially Expressed Genes", subtitle="Medial Amygdala") + 
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
ggsave(plot=manhattan_MeA,file="png/Manhattan_MedialAmyg_PTSD_labels_adjP.png",width=16,height=12,dpi=320)


#####################################
#dACC
#####################################

#without minor axis
manhattan_dACC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_dACC, x=building$gene_loc, color=building$dACC_color),size=3) +
scale_colour_manual(values=dacc_col_values) +
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


#####################################
#DLPFC
#####################################

#without minor axis
manhattan_DLPFC <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_DLPFC, x=building$gene_loc, color=building$DLPFC_color),size=3) +
scale_colour_manual(values=dlpfc_col_values) +
scale_y_continuous(limits=c(-1.5,1),labels = c(seq(1.5, 0, by = -0.25),seq(0.25, 1, by = 0.25)), breaks=seq(-1.5, 1, by = 0.25)) +
geom_hline(yintercept=1, color="#B23AEE",size=0.75) +
geom_hline(yintercept=-1,color="#B23AEE",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = subset(building, adj.P.Val_PTSD_DLPFC <0.1), aes(label = Symbol, y=final.P.adj_DLPFC, x=gene_loc),color=text_color, 
	family=font,fontface='bold',size=6,min.segment.length = unit(0, 'lines'),nudge_x=50000000,hjust=0,vjust=0.5) +
ggtitle("Differentially Expressed Genes", subtitle="Dorsolateral Prefrontal Cortex") + 
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
ggsave(plot=manhattan_DLPFC,file="png/Manhattan_DLPFC_PTSD_labels_adjP.png",width=16,height=12,dpi=320)


#####################################
#Pan
#####################################

#without minor axis
manhattan_Pan <- ggplot() +
geom_point(data=building,aes(y=building$final.P.adj_Pan, x=building$gene_loc, color=building$Pan_color),size=3) +
scale_colour_manual(values=pan_col_values) +
scale_y_continuous(limits=c(-3.75,3.25),labels = c(seq(3.5, 0, by = -0.5),seq(0.5, 3.5, by = 0.5)), breaks=seq(-3.5, 3.5, by = 0.5)) +
geom_hline(yintercept=1, color="#CD2626",size=0.75) +
geom_hline(yintercept=-1,color="#CD2626",size=0.75) +
geom_hline(yintercept=0, color="#222222",linetype="dashed",size=1,alpha=0.5) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][c(1,4,5,10),], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][7,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,hjust=0,vjust=0,nudge_y=0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][2,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][3,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,hjust=0,vjust=0,nudge_y=0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][6,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,hjust=0,vjust=0,nudge_y=0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][8,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,hjust=0,vjust=1,nudge_y=-0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan,decreasing=TRUE),][9,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_x = -55000000,hjust=1,vjust=1) +
geom_text_repel(data = building[order(building$final.P.adj_Pan),][c(2:6,8,10),], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.3) +
geom_text_repel(data = building[order(building$final.P.adj_Pan),][9,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y=-0.3,hjust=0,vjust=1) +
geom_text_repel(data = building[order(building$final.P.adj_Pan),][1,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,hjust=0,vjust=0,nudge_y=0.25) +
geom_text_repel(data = building[order(building$final.P.adj_Pan),][7,], aes(label = Symbol, y=final.P.adj_Pan, x=gene_loc),min.segment.length = unit(0, 'lines'),color=text_color, 
	family=font,fontface='bold',size=6,nudge_y = 0.45) +
ggtitle("Differentially Expressed Genes", subtitle="Four Region Interaction Model, PTSD Term") + 
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
ggsave(plot=manhattan_Pan,file="png/Manhattan_4RegionIntxn_GroupPTSD_PTSD_labels_adjP.png",width=16,height=12,dpi=320)
