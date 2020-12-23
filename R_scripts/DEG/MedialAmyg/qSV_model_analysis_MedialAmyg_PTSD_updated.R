#####################################
#Make qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/make_qSVs.R)
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)

#load rse_gene object
load('../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("../rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])


## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Find genetic race principal components that associate with Group to use in the model below
summary(lm(rse_gene$snpPC1 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC2 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC3 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC8 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC9 ~ rse_gene$Group))
summary(lm(rse_gene$snpPC10 ~ rse_gene$Group))
#These (above) were the components significantly associated with one or both groups (MDD or PTSD)

## model
with(colData(rse_gene), table(Group))

#> with(colData(rse_gene), table(Group))
#Group
#Control     MDD    PTSD 
#    431     432     422 

## also add snpPCs based on association with group
mod = model.matrix(~Group + AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

#load cov_rse object from PTSD data
load("../rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

#Make sure samples line up across qSVs and gene counts
cov_rse = cov_rse[,colnames(rse_gene)]
#Check
identical(colnames(cov_rse), colnames(rse_gene))

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod) # 19
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

#> getPcaVars(qsvBonf)[1:k]
#[1] 75.500  5.300  3.170  1.760  1.230  1.020  0.816  0.407  0.377  0.352
#[11]  0.326  0.290  0.269  0.242  0.217  0.207  0.189  0.161  0.154

modQsva = cbind(mod, qSVs)

#save(qsvBonf, qSVs, mod, modQsva, file = '../rdas/PTSD_qsvs.Rdata')
#dir.create("pdf", showWarnings = FALSE)
#pdf(file='../pdf/qsvs_var_explained.pdf', useDingbats = FALSE)
#plot(getPcaVars(qsvBonf)[1:k], pch=20)
#dev.off()


#####################################
#Explore qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_qsvs.R)
#####################################

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)
library('readxl')
library('devtools')

# gene-level expression for degradation data
#If running this with above code, then don't need to reload objects


#load('../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
## add MDS (get ethnicity via genotype)
#load("../rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
#rownames(mds) = ss(rownames(mds),"_")
#colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])
#rse_gene <- rse_gene[, rse_gene$Region %in% c('DLPFC', 'BasoAmyg', 'MedialAmyg', 'dACC')]
#gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
#rse_gene <- rse_gene[gIndex , ]

getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))

# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
# > getPcaVars(pca)[1:5]
#20.40 11.80  6.78  5.13  2.69


## degradation plots
pdf('../pdf/degradation.pdf', useDingbats = FALSE)
plot(pca$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate", pch=20,
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,2] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pc2: ",getPcaVars(pca)[2],"% Var Expl"))
plot(pca$x[,1] ~ pd$RIN,
     xlab = "RIN", pch=20,
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
   ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$overallMapRate,
     xlab = "Overall Map Rate", pch=20,
   ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()



#PTSD plots
pd$Group = factor(pd$Group)
pd$Sex = factor(pd$Sex)
pd$Race = factor(pd$Race)

pdf('../pdf/PTSD_plots.pdf', useDingbats = FALSE)

plot(pca$x[,1] ~ pd$Group,
     xlab = "Group",
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$Sex,
     xlab = "Sex",
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$AgeDeath,
     xlab = "Age", pch=20,
     ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')

plot(pca$x[,2] ~ pd$AgeDeath,
     xlab = "Age", pch=20,
     ylab=paste0("pc2: ",getPcaVars(pca)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')

plot(pca$x[,1] ~ pd$Race,
     xlab = "Race", pch=20,
	 ylab=paste0("pc1: ",getPcaVars(pca)[1],"% Var Expl"))


plot(pca$x[,2] ~ pd$Race,
     xlab = "Race", pch=20,
	 ylab=paste0("pc2: ",getPcaVars(pca)[2],"% Var Expl"))


dev.off()


# top 1000 expressed regions associated with degradation
#load("../rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

## get qSVs for top bonferroni
#qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
#getPcaVars(qsvBonf)[1:5]
# > getPcaVars(qsvBonf)[1:5]
#[1] 75.50  5.30  3.17  1.76  1.23

##PTSD qsv plots
pdf("../pdf/qSVs_PTSD.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$Group,
     xlab = "Group",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$AgeDeath,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,2] ~ pd$AgeDeath,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()


#####################################
#Case-Control Medial Amygdala (https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC.R)
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')
library(recount)

#If following from above, don't need to reload objects

#load rse object from PTSD data
#load('../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata', verbose = TRUE)

## add MDS (get ethnicity via genotype)
#load("../rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
#rownames(mds) = ss(rownames(mds),"_")
#colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff
##		do across all regions so we're looking at the same genes
#gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
#rse_gene <- rse_gene[gIndex , ]

#load qsvBonf, qSVs, mod, modQsva object from PTSD data
#load("../rdas/PTSD_qsvs.Rdata", verbose = TRUE)


###Analysis

#filter for MedialAmyg
keepIndex = which(rse_gene$Region == "MedialAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <- !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[c(1, 18:36)] = c("Int", gsub("PC", "qSV", colnames(modQsva)[18:36]))

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('../pdf/medialamyg_voom_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()

fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#MDD vs controls (analysis 1)
sigGeneMDD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneMDD) = paste0(colnames(sigGeneMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigGenePTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGenePTSD) = paste0(colnames(sigGenePTSD), "_PTSD")

#Interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")

## no qSVA
pdf('../pdf/medialamyg_voom_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)

sigGene0MDD = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene0MDD) = paste0(colnames(sigGene0MDD), "_0MDD")

sigGene0PTSD = topTable(eBGene0,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene0PTSD) = paste0(colnames(sigGene0PTSD), "_0PTSD")

sigGene0Dx = topTable(eBGene0,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene0Dx) = paste0(colnames(sigGene0Dx), "_0ANOVA")

## no adjustment vars
pdf('../pdf/medialamyg_voom_noadj.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Group)), plot=TRUE)
dev.off()

fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)

sigGeneNoAdjMDD = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneNoAdjMDD) = paste0(colnames(sigGeneNoAdjMDD), "_NoAdjMDD")

sigGeneNoAdjPTSD = topTable(eBGeneNoAdj,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneNoAdjPTSD) = paste0(colnames(sigGeneNoAdjPTSD), "_NoAdjPTSD")

sigGeneNoAdjDx = topTable(eBGeneNoAdj,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneNoAdjDx) = paste0(colnames(sigGeneNoAdjDx), "_NoAdjANOVA")


#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse object from PTSD data
load('../rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata', verbose = TRUE)

## add MDS (get ethnicity via genotype)
load("../rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

#Redefine model using PTSD only

mod = model.matrix(~onlyPTSD+ AgeDeath + Sex + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

#load cov_rse object from PTSD data
load("../rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

#Make sure samples line up across qSVs and gene counts
cov_rse = cov_rse[,colnames(rse_gene)]
#Check
identical(colnames(cov_rse), colnames(rse_gene))

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod) 
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

#> getPcaVars(qsvBonf)[1:k] # 19
#[1] 75.500  5.300  3.170  1.760  1.230  1.020  0.816  0.407  0.377  0.352
#[11]  0.326  0.290  0.269  0.242  0.217  0.207  0.189  0.161  0.154 

modQsva = cbind(mod, qSVs)

#save(qsvBonf, qSVs, mod, modQsva, file = '../rdas/PTSD_onlyPTSD_qsvs.Rdata')
#pdf(file='../pdf/qsvs_var_explained_onlyPTSD.pdf', useDingbats = FALSE)
#plot(getPcaVars(qsvBonf)[1:k], pch=20)
#dev.off()

#filter for MedialAmyg
keepIndex = which(rse_gene$Region == "MedialAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <- !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('../pdf/medialamyg_voom_onlyPTSD_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()

fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")

## no qSVA
pdf('../pdf/medialamyg_voom_onlyPTSD_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()

fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)

sigGene0onlyPTSD = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene0onlyPTSD) = paste0(colnames(sigGene0onlyPTSD), "_0onlyPTSD")

## no adjustment vars
pdf('../pdf/medialamyg_voom_onlyPTSD_noadj.pdf', useDingbats = FALSE)
vGeneNoAdjonlyPTSD = voom(dge, with(colData(rse_gene), model.matrix( ~ onlyPTSD)), plot=TRUE)
dev.off()

fitGeneNoAdjonlyPTSD = lmFit(vGeneNoAdjonlyPTSD)
eBGeneNoAdjonlyPTSD = eBayes(fitGeneNoAdjonlyPTSD)

sigGeneNoAdjonlyPTSD = topTable(eBGeneNoAdjonlyPTSD,coef=2,
	p.value = 1,number=nrow(rse_gene))
colnames(sigGeneNoAdjonlyPTSD) = paste0(colnames(sigGeneNoAdjonlyPTSD), "_NoAdjonlyPTSD")


###Stuff we want
#Merge qSVA into one big genestats dataframe
geneStats = cbind(sigGenePTSD[,c(11,13:15)], sigGeneMDD[,c(11,13:15)],sigGeneDx[,c(14:16)],sigGeneonlyPTSD[,c(11,13:15)], PTSDvsMDDPost[,c(11,13:15)])
geneStats = cbind(geneStats, rowData(rse_gene))

## write out qSVA-based stats
geneStats_medialamyg = geneStats
save(geneStats_medialamyg, file = "../rdas/geneStats_DE_qSVA_MedialAmyg_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_medialamyg, file = gzfile("../csvs/geneStats_DE_qSVA_MedialAmyg_threeGroup.csv.gz"))

#Merge noqSVA into one big genestats dataframe
geneStats = cbind(sigGene0PTSD[,c(11,13:15)], sigGene0MDD[,c(11,13:15)],sigGene0Dx[,c(14:16)],sigGene0onlyPTSD[,c(11,13:15)])
geneStats = cbind(geneStats, rowData(rse_gene))

## write out no qSVA-based stats
geneStats_medialamyg = geneStats
save(geneStats_medialamyg, file = "../rdas/geneStats_DE_noqSVA_MedialAmyg_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_medialamyg, file = gzfile("../csvs/geneStats_DE_noqSVA_MedialAmyg_threeGroup.csv.gz"))

#Merge no adjustment vars into one big genestats dataframe
geneStats = cbind(sigGeneNoAdjPTSD[,c(11,13:15)], sigGeneNoAdjMDD[,c(11,13:15)],sigGeneNoAdjDx[,c(14:16)],sigGeneNoAdjonlyPTSD[,c(11,13:15)])
geneStats = cbind(geneStats, rowData(rse_gene))

## write out no adjustment vars-based stats
geneStats_medialamyg = geneStats
save(geneStats_medialamyg, file = "../rdas/geneStats_DE_NoAdj_MedialAmyg_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_medialamyg, file = gzfile("../csvs/geneStats_DE_NoAdj_MedialAmyg_threeGroup.csv.gz"))




#####################################
#Explore Case-Control MedialAmyg (https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_case_control.R)
#####################################






## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

#> print('Reproducibility information:')
#[1] "Reproducibility information:"
#> Sys.time()
#[1] "2019-01-30 10:54:14 EST"
#> proc.time()
#    user   system  elapsed 
# 633.510    5.313 1060.777 
#> options(width = 120)
#> session_info()
#─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                                 
# version  R version 3.5.1 Patched (2018-10-29 r75535)           
# os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
# system   x86_64, linux-gnu                                     
# ui       X11                                                   
# language (EN)                                                  
# collate  en_US.UTF-8                                           
# ctype    en_US.UTF-8                                           
# tz       US/Eastern                                            
# date     2019-01-30                                            

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 package              * version   date       lib source                                   
 acepack                1.4.1     2016-10-29 [2] CRAN (R 3.5.0)                           
 annotate               1.60.0    2018-10-30 [2] Bioconductor                             
 AnnotationDbi          1.44.0    2018-10-30 [2] Bioconductor                             
 assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)                           
 backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)                           
 base64enc              0.1-3     2015-07-28 [2] CRAN (R 3.5.0)                           
 bibtex                 0.4.2     2017-06-30 [2] CRAN (R 3.5.1)                           
 bindr                  0.1.1     2018-03-13 [2] CRAN (R 3.5.0)                           
 bindrcpp               0.2.2     2018-03-29 [2] CRAN (R 3.5.0)                           
 Biobase              * 2.42.0    2018-10-30 [2] Bioconductor                             
 BiocGenerics         * 0.28.0    2018-10-30 [2] Bioconductor                             
 BiocParallel         * 1.16.5    2019-01-04 [2] Bioconductor                             
 biomaRt                2.38.0    2018-10-30 [2] Bioconductor                             
 Biostrings             2.50.1    2018-11-06 [2] Bioconductor                             
 bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.0)                           
 bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)                           
 bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)                           
 blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)                           
 BSgenome               1.50.0    2018-10-30 [2] Bioconductor                             
 bumphunter             1.24.5    2018-12-01 [2] Bioconductor                             
 callr                  3.1.1     2018-12-21 [2] CRAN (R 3.5.1)                           
 cellranger             1.1.0     2016-07-27 [2] CRAN (R 3.5.0)                           
 checkmate              1.9.1     2019-01-15 [2] CRAN (R 3.5.1)                           
 cli                    1.0.1     2018-09-25 [2] CRAN (R 3.5.1)                           
 cluster                2.0.7-1   2018-04-13 [3] CRAN (R 3.5.1)                           
 codetools              0.2-15    2016-10-05 [3] CRAN (R 3.5.1)                           
 colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)                           
 crayon                 1.3.4     2017-09-16 [2] CRAN (R 3.5.0)                           
 data.table             1.12.0    2019-01-13 [2] CRAN (R 3.5.1)                           
 DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)                           
 DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor                             
 derfinder              1.16.1    2018-12-03 [2] Bioconductor                             
 derfinderHelper        1.16.1    2018-12-03 [2] Bioconductor                             
 desc                   1.2.0     2018-05-01 [2] CRAN (R 3.5.1)                           
 devtools             * 2.0.1     2018-10-26 [2] CRAN (R 3.5.1)                           
 digest                 0.6.18    2018-10-10 [2] CRAN (R 3.5.1)                           
 doRNG                  1.7.1     2018-06-22 [2] CRAN (R 3.5.1)                           
 downloader             0.4       2015-07-09 [2] CRAN (R 3.5.0)                           
 dplyr                  0.7.8     2018-11-10 [2] CRAN (R 3.5.1)                           
 edgeR                * 3.24.3    2019-01-02 [2] Bioconductor                             
 foreach                1.4.4     2017-12-12 [2] CRAN (R 3.5.0)                           
 foreign                0.8-71    2018-07-20 [3] CRAN (R 3.5.1)                           
 Formula                1.2-3     2018-05-03 [2] CRAN (R 3.5.0)                           
 fs                     1.2.6     2018-08-23 [2] CRAN (R 3.5.1)                           
 genefilter           * 1.64.0    2018-10-30 [2] Bioconductor                             
 GenomeInfoDb         * 1.18.1    2018-11-12 [2] Bioconductor                             
 GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor                             
 GenomicAlignments      1.18.1    2019-01-04 [2] Bioconductor                             
 GenomicFeatures        1.34.1    2018-11-03 [2] Bioconductor                             
 GenomicFiles           1.18.0    2018-10-30 [2] Bioconductor                             
 GenomicRanges        * 1.34.0    2018-10-30 [2] Bioconductor                             
 GEOquery               2.50.4    2018-12-13 [2] Bioconductor                             
 ggplot2                3.1.0     2018-10-25 [2] CRAN (R 3.5.1)                           
 glue                   1.3.0     2018-07-17 [2] CRAN (R 3.5.1)                           
 gridExtra              2.3       2017-09-09 [2] CRAN (R 3.5.0)                           
 gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)                           
 Hmisc                  4.1-1     2018-01-03 [2] CRAN (R 3.5.0)                           
 hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)                           
 htmlTable              1.13.1    2019-01-07 [2] CRAN (R 3.5.1)                           
 htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)                           
 htmlwidgets            1.3       2018-09-30 [2] CRAN (R 3.5.1)                           
 httr                   1.4.0     2018-12-11 [2] CRAN (R 3.5.1)                           
 IRanges              * 2.16.0    2018-10-30 [2] Bioconductor                             
 iterators              1.0.10    2018-07-13 [2] CRAN (R 3.5.1)                           
 jaffelab             * 0.99.22   2018-07-23 [1] Github (LieberInstitute/jaffelab@a9e7377)
 jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)                           
 knitr                  1.21      2018-12-10 [2] CRAN (R 3.5.1)                           
 lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)                           
 latticeExtra           0.6-28    2016-02-09 [2] CRAN (R 3.5.0)                           
 lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)                           
 limma                * 3.38.3    2018-12-02 [2] Bioconductor                             
 locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)                           
 magrittr               1.5       2014-11-22 [2] CRAN (R 3.5.0)                           
 Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)                           
 matrixStats          * 0.54.0    2018-07-23 [2] CRAN (R 3.5.1)                           
 memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)                           
 mgcv                 * 1.8-26    2018-11-21 [3] CRAN (R 3.5.1)                           
 munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)                           
 nlme                 * 3.1-137   2018-04-07 [3] CRAN (R 3.5.1)                           
 nnet                   7.3-12    2016-02-02 [3] CRAN (R 3.5.1)                           
 pillar                 1.3.1     2018-12-15 [2] CRAN (R 3.5.1)                           
 pkgbuild               1.0.2     2018-10-16 [2] CRAN (R 3.5.1)                           
 pkgconfig              2.0.2     2018-08-16 [2] CRAN (R 3.5.1)                           
 pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.5.1)                           
 pkgmaker               0.27      2018-05-25 [2] CRAN (R 3.5.0)                           
 plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)                           
 prettyunits            1.0.2     2015-07-13 [2] CRAN (R 3.5.0)                           
 processx               3.2.1     2018-12-05 [2] CRAN (R 3.5.1)                           
 progress               1.2.0     2018-06-14 [2] CRAN (R 3.5.1)                           
 ps                     1.3.0     2018-12-21 [2] CRAN (R 3.5.1)                           
 purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)                           
 qvalue                 2.14.1    2019-01-10 [2] Bioconductor                             
 R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)                           
 rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)                           
 RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)                           
 Rcpp                   1.0.0     2018-11-07 [2] CRAN (R 3.5.1)                           
 RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)                           
 readr                  1.3.0     2018-12-11 [2] CRAN (R 3.5.1)                           
 readxl               * 1.2.0     2018-12-19 [2] CRAN (R 3.5.1)                           
 recount              * 1.8.1     2018-12-03 [1] Bioconductor                             
 registry               0.5       2017-12-03 [2] CRAN (R 3.5.0)                           
 remotes                2.0.2     2018-10-30 [2] CRAN (R 3.5.1)                           
 rentrez                1.2.1     2018-03-05 [2] CRAN (R 3.5.0)                           
 reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)                           
 rlang                  0.3.1     2019-01-08 [2] CRAN (R 3.5.1)                           
 rngtools               1.3.1     2018-05-15 [2] CRAN (R 3.5.0)                           
 rpart                  4.1-13    2018-02-23 [3] CRAN (R 3.5.1)                           
 rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.5.0)                           
 Rsamtools              1.34.0    2018-10-30 [2] Bioconductor                             
 RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)                           
 rstudioapi             0.9.0     2019-01-09 [2] CRAN (R 3.5.1)                           
 rtracklayer            1.42.1    2018-11-22 [2] Bioconductor                             
 S4Vectors            * 0.20.1    2018-11-09 [2] Bioconductor                             
 scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)                           
 segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)                           
 sessioninfo            1.1.1     2018-11-05 [2] CRAN (R 3.5.1)                           
 stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)                           
 stringr                1.3.1     2018-05-10 [2] CRAN (R 3.5.1)                           
 SummarizedExperiment * 1.12.0    2018-10-30 [2] Bioconductor                             
 survival               2.43-3    2018-11-26 [3] CRAN (R 3.5.1)                           
 sva                  * 3.30.1    2019-01-04 [2] Bioconductor                             
 tibble                 2.0.0     2019-01-04 [2] CRAN (R 3.5.1)                           
 tidyr                  0.8.2     2018-10-28 [2] CRAN (R 3.5.1)                           
 tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)                           
 usethis              * 1.4.0     2018-08-14 [2] CRAN (R 3.5.1)                           
 VariantAnnotation      1.28.8    2019-01-10 [2] Bioconductor                             
 withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)                           
 xfun                   0.4       2018-10-23 [2] CRAN (R 3.5.1)                           
 XML                    3.98-1.16 2018-08-19 [2] CRAN (R 3.5.1)                           
 xml2                   1.2.0     2018-01-24 [2] CRAN (R 3.5.0)                           
 xtable                 1.8-3     2018-08-29 [2] CRAN (R 3.5.1)                           
 XVector                0.22.0    2018-10-30 [2] Bioconductor                             
 zlibbioc               1.28.0    2018-10-30 [2] Bioconductor                             

[1] /users/bbarry/R/x86_64-pc-linux-gnu-library/3.5
[2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
[3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library


