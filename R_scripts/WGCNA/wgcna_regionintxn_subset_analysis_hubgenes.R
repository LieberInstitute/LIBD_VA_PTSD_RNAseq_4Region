#Trying my hand at determining hub genes with https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#especially https://pdfs.semanticscholar.org/5e42/e2185c54874277794395e5825808e5f5709c.pdf as my main reference
#key function: signedKME()
#also using Peter Langfelder's replies here for guidance: https://support.bioconductor.org/p/69858/


##########################################################################################################
#####Amyg
##########################################################################################################

################
#PTSD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(devtools)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_PTSD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Amyg_PTSD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_PTSD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,3,4,5,6,20,2,7:19,21:40)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
##Define hub genes as having high gene significance (association of gene with disease) and high intramodular connectivity 
	#(define a module eigengene-based connectivity measure for each gene as the correlation between a the gene expression and the module eigengene.)

#####High gene significance:

#We've already completed our differential expression analysis, identifying genes DE between cases and controls (ie, high gene significance)
#because we have ~2 brain regions from each donor, we'd want to incorporate duplicate correlation so i wouldnt just use cor like in the examples in the tutorial
load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_PTSD <- geneStatsall[geneStatsall$P.Value_PTSD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,]
#modules 2 (GO = immune/inflam terms), 18(GO = extracell matrix/cell surface proteins, axon guidance), 20 (GO = synaptic vesicle membrane)

sig_genes_connect2 <- rownames(datKME)[abs(datKME$MM.2) > 0.5]
sig_genes_connect17 <- rownames(datKME)[abs(datKME$MM.17) > 0.5]
sig_genes_connect18 <- rownames(datKME)[abs(datKME$MM.18) > 0.5]
sig_genes_connect20 <- rownames(datKME)[abs(datKME$MM.20) > 0.5]


intx2 <- intersect(sig_genes_PTSD, sig_genes_connect2)
intx17 <- intersect(sig_genes_PTSD, sig_genes_connect17)
intx18 <- intersect(sig_genes_PTSD, sig_genes_connect18)
intx20 <- intersect(sig_genes_PTSD, sig_genes_connect20)

overlap<-datKME[rownames(datKME) %in% sig_genes_PTSD,]

save(coefAdj,sig_modules,sig_genes_connect2,sig_genes_connect17,sig_genes_connect18,sig_genes_connect20,
	intx2,intx17,intx18,intx20,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Amyg_PTSD_hubgenes.rda")

################
#MDD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_MDD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Amyg_MDD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_MDD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

## model
load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,2,4,5,6,19,3,7:18,20:40)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
#####High gene significance:

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_MDD <- geneStatsall[geneStatsall$P.Value_MDD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,]

sig_genes_connect1 <- rownames(datKME)[abs(datKME$MM.1) > 0.5]
sig_genes_connect10 <- rownames(datKME)[abs(datKME$MM.10) > 0.5]
sig_genes_connect20 <- rownames(datKME)[abs(datKME$MM.20) > 0.5]
sig_genes_connect22 <- rownames(datKME)[abs(datKME$MM.22) > 0.5]

intx1 <- intersect(sig_genes_MDD, sig_genes_connect1)
intx10 <- intersect(sig_genes_MDD, sig_genes_connect10)
intx20 <- intersect(sig_genes_MDD, sig_genes_connect20)
intx22 <- intersect(sig_genes_MDD, sig_genes_connect22)

overlap<-datKME[rownames(datKME) %in% sig_genes_MDD,]

save(coefAdj,sig_modules,sig_genes_connect1,sig_genes_connect10,sig_genes_connect20,sig_genes_connect22,
	intx1,intx10,intx20,intx22,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Amyg_MDD_hubgenes.rda")

################
#onlyPTSD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_onlyPTSD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Amyg_onlyPTSD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_onlyPTSD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

## model #

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD_Amyg.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1:5,18,6:17,19:38)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
#####High gene significance:

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_onlyPTSD <- geneStatsall[geneStatsall$P.Value_onlyPTSD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,,drop=FALSE]

#modules 2 (GO = immune response/inflam)

sig_genes_connect2 <- rownames(datKME)[abs(datKME$MM.2) > 0.5]

intx2 <- intersect(sig_genes_onlyPTSD, sig_genes_connect2)

overlap<-datKME[rownames(datKME) %in% sig_genes_onlyPTSD,]

save(coefAdj,sig_modules,sig_genes_connect2,intx2,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Amyg_onlyPTSD_hubgenes.rda")


##########################################################################################################
#####Cortex
##########################################################################################################

################
#PTSD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_PTSD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Cortex_PTSD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_PTSD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Cortex")
rse_gene <- rse_gene[, gkeepIndex]

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
load('rdas/PTSD_qsvs_Regionintxn_Cortex.Rdata',verbose=TRUE)
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,3,4,5,6,20,2,7:19,21:40)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
#####High gene significance:
load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_PTSD <- geneStatsall[geneStatsall$P.Value_PTSD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,]

sig_genes_connect12 <- rownames(datKME)[abs(datKME$MM.12) > 0.5]
sig_genes_connect1 <- rownames(datKME)[abs(datKME$MM.1) > 0.5]
sig_genes_connect15 <- rownames(datKME)[abs(datKME$MM.15) > 0.5]
sig_genes_connect10 <- rownames(datKME)[abs(datKME$MM.10) > 0.5]
sig_genes_connect11 <- rownames(datKME)[abs(datKME$MM.11) > 0.5]
sig_genes_connect17 <- rownames(datKME)[abs(datKME$MM.17) > 0.5]
sig_genes_connect3 <- rownames(datKME)[abs(datKME$MM.3) > 0.5]
sig_genes_connect7 <- rownames(datKME)[abs(datKME$MM.7) > 0.5]
sig_genes_connect8 <- rownames(datKME)[abs(datKME$MM.8) > 0.5]
sig_genes_connect2 <- rownames(datKME)[abs(datKME$MM.2) > 0.5]
sig_genes_connect0 <- rownames(datKME)[abs(datKME$MM.0) > 0.5]

intx12 <- intersect(sig_genes_PTSD, sig_genes_connect12)
intx1 <- intersect(sig_genes_PTSD, sig_genes_connect1)
intx15 <- intersect(sig_genes_PTSD, sig_genes_connect15)
intx10 <- intersect(sig_genes_PTSD, sig_genes_connect10)
intx11 <- intersect(sig_genes_PTSD, sig_genes_connect11)
intx17 <- intersect(sig_genes_PTSD, sig_genes_connect17)
intx3 <- intersect(sig_genes_PTSD, sig_genes_connect3)
intx7 <- intersect(sig_genes_PTSD, sig_genes_connect7)
intx8 <- intersect(sig_genes_PTSD, sig_genes_connect8)
intx2 <- intersect(sig_genes_PTSD, sig_genes_connect2)
intx0 <- intersect(sig_genes_PTSD, sig_genes_connect0)


overlap<-datKME[rownames(datKME) %in% sig_genes_PTSD,]

save(coefAdj,sig_modules,sig_genes_connect12,sig_genes_connect1,sig_genes_connect15,sig_genes_connect10,sig_genes_connect11,sig_genes_connect17,sig_genes_connect3,sig_genes_connect7,
	sig_genes_connect8,sig_genes_connect2,sig_genes_connect0,intx12,intx1,intx15,intx10,intx11,intx17,intx3,intx7,intx8,intx2,intx0
	,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Cortex_PTSD_hubgenes.rda")

################
#MDD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_MDD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Cortex_MDD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_MDD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Cortex")
rse_gene <- rse_gene[, gkeepIndex]

## model
load('rdas/PTSD_qsvs_Regionintxn_Cortex.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1,2,4,5,6,19,3,7:18,20:40)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
#####High gene significance:

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_MDD <- geneStatsall[geneStatsall$P.Value_MDD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,]

sig_genes_connect2<-rownames(datKME)[abs(datKME$MM.2) > 0.5]  
sig_genes_connect9<-rownames(datKME)[abs(datKME$MM.9) > 0.5]  
sig_genes_connect5<-rownames(datKME)[abs(datKME$MM.5) > 0.5]  
sig_genes_connect17<-rownames(datKME)[abs(datKME$MM.17) > 0.5]
sig_genes_connect20<-rownames(datKME)[abs(datKME$MM.20) > 0.5]
sig_genes_connect7<-rownames(datKME)[abs(datKME$MM.7) > 0.5]  
sig_genes_connect10<-rownames(datKME)[abs(datKME$MM.10) > 0.5]
sig_genes_connect27<-rownames(datKME)[abs(datKME$MM.27) > 0.5]
sig_genes_connect6<-rownames(datKME)[abs(datKME$MM.6) > 0.5]  
sig_genes_connect8<-rownames(datKME)[abs(datKME$MM.8) > 0.5]  
sig_genes_connect21<-rownames(datKME)[abs(datKME$MM.21) > 0.5]
sig_genes_connect4<-rownames(datKME)[abs(datKME$MM.4) > 0.5]  
sig_genes_connect11<-rownames(datKME)[abs(datKME$MM.11) > 0.5]
sig_genes_connect12<-rownames(datKME)[abs(datKME$MM.12) > 0.5]
sig_genes_connect22<-rownames(datKME)[abs(datKME$MM.22) > 0.5]
sig_genes_connect0<-rownames(datKME)[abs(datKME$MM.0) > 0.5]

intx2<-intersect(sig_genes_MDD, sig_genes_connect2)  
intx9<-intersect(sig_genes_MDD, sig_genes_connect9)  
intx5<-intersect(sig_genes_MDD, sig_genes_connect5)  
intx17<-intersect(sig_genes_MDD, sig_genes_connect17)
intx20<-intersect(sig_genes_MDD, sig_genes_connect20)
intx7<-intersect(sig_genes_MDD, sig_genes_connect7)  
intx10<-intersect(sig_genes_MDD, sig_genes_connect10)
intx27<-intersect(sig_genes_MDD, sig_genes_connect27)
intx6<-intersect(sig_genes_MDD, sig_genes_connect6)  
intx8<-intersect(sig_genes_MDD, sig_genes_connect8)  
intx21<-intersect(sig_genes_MDD, sig_genes_connect21)
intx4<-intersect(sig_genes_MDD, sig_genes_connect4)  
intx11<-intersect(sig_genes_MDD, sig_genes_connect11)
intx12<-intersect(sig_genes_MDD, sig_genes_connect12)
intx22<-intersect(sig_genes_MDD, sig_genes_connect22)
intx0<-intersect(sig_genes_MDD, sig_genes_connect0)

overlap<-datKME[rownames(datKME) %in% sig_genes_MDD,]

save(coefAdj,sig_modules, sig_genes_connect2, sig_genes_connect9, sig_genes_connect5,sig_genes_connect17, sig_genes_connect20, sig_genes_connect7,
	sig_genes_connect10, sig_genes_connect27, sig_genes_connect6, sig_genes_connect8,  sig_genes_connect21,sig_genes_connect4, 
	sig_genes_connect11, sig_genes_connect12, sig_genes_connect22,sig_genes_connect0, intx2, intx9, intx5, intx17, intx20, intx7, intx10, intx27, intx6, 
	intx8, intx21, intx4, intx11, intx12, intx22, intx0,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Cortex_MDD_hubgenes.rda")


################
#onlyPTSD
################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(clusterProfiler)
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#load objects
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_onlyPTSD.rda",verbose=TRUE)
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
load("rdas/WGCNA/MEvsDx_Cortex_onlyPTSD.rda",verbose=TRUE)
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_onlyPTSD.rda",verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

# filter ##
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")
gkeepIndex = which(rse_gene$CatRegion == "Cortex")
rse_gene <- rse_gene[, gkeepIndex]

## model #

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD_Cortex.Rdata',verbose=TRUE)

##ensure intercept, region, interaction, group, agedeath, and sex are up first for cleaning
colnames(modQsva)
dim(modQsva)
modQsva = modQsva[,c(1:5,18,6:17,19:38)]
head(modQsva)

## clean expression
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
geneExprsClean = cleaningY(geneExprs, modQsva, P=6)

#####Finding genes with high gene significance and high intramodular connectivity in interesting modules
#####High gene significance:

load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",verbose=TRUE)
identical(rownames(geneExprsClean),geneStatsall$gencodeID,attrib.as.set=FALSE)

sig_genes_onlyPTSD <- geneStatsall[geneStatsall$P.Value_onlyPTSD < 0.005,"gencodeID"]

#####High intramodular connectivity (note, the gene doesnt have to be a member of that module)
datME <- net$MEs #Rows correspond to samples and columns to module eigengenes.

datExpr <- t(geneExprsClean) #Rows correspond to samples and columns to genes.

datKME=signedKME(datExpr, datME, outputColumnName="MM.")

#####Interesting modules
sig_modules<-coefAdj[coefAdj[,5] < 0.05,,drop=FALSE]


sig_genes_connect1 <- rownames(datKME)[abs(datKME$MM.1) > 0.5]
sig_genes_connect3 <- rownames(datKME)[abs(datKME$MM.3) > 0.5]

intx1 <- intersect(sig_genes_onlyPTSD, sig_genes_connect1)
intx3 <- intersect(sig_genes_onlyPTSD, sig_genes_connect3)

overlap<-datKME[rownames(datKME) %in% sig_genes_onlyPTSD,]

save(coefAdj,sig_modules,sig_genes_connect1,sig_genes_connect3,intx1,intx3,overlap,file="rdas/WGCNA/wgcna_regionintxn_subset_Cortex_onlyPTSD_hubgenes.rda")

###############################################################################################################################################
