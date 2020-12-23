#####Stratify by sex analysis#######


############# Males ##############

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]


load('rdas/PTSD_qsvs.Rdata')


#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by sex
malekeepIndex = which(rse_gene$Sex == "M")
rse_gene <- rse_gene[, malekeepIndex]

#Get rid of sex columns
colIndex <-  !grepl("Sex", colnames(mod))
mod <- mod[malekeepIndex, colIndex]
modQsva <- modQsva[malekeepIndex,!grepl("Sex", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
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
PTSDvsMDDPostGene = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPostGene) = paste0(colnames(PTSDvsMDDPostGene), "_PTSDvsMDD")




#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs.Rdata')

#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by sex
malekeepIndex = which(rse_gene$Sex == "M")
rse_gene <- rse_gene[, malekeepIndex]

#Get rid of sex columns
colIndex <-  !grepl("Sex", colnames(mod))
mod <- mod[malekeepIndex, colIndex]
modQsva <- modQsva[malekeepIndex,!grepl("Sex", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


###Stuff we want

###### GENE ########


#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPostGene)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
geneStats_BasoAmygall_male = geneStatsall
save(geneStats_BasoAmygall_male, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_males_threeGroup.rda")

write.csv(geneStats_BasoAmygall_male, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_males_threeGroup.csv.gz"))





############# Females ##############



library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]


load('rdas/PTSD_qsvs.Rdata')


#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by sex
femalekeepIndex = which(rse_gene$Sex == "F")
rse_gene <- rse_gene[, femalekeepIndex]

#Get rid of sex columns
colIndex <-  !grepl("Sex", colnames(mod))
mod <- mod[femalekeepIndex, colIndex]
modQsva <- modQsva[femalekeepIndex,!grepl("Sex", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
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
PTSDvsMDDPostGene = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPostGene) = paste0(colnames(PTSDvsMDDPostGene), "_PTSDvsMDD")




#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')


## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs.Rdata')

#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


#stratify by sex
femalekeepIndex = which(rse_gene$Sex == "F")
rse_gene <- rse_gene[, femalekeepIndex]

#Get rid of sex columns
colIndex <-  !grepl("Sex", colnames(mod))
mod <- mod[femalekeepIndex, colIndex]
modQsva <- modQsva[femalekeepIndex,!grepl("Sex", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


###Stuff we want

###### GENE ########


#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPostGene)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
geneStats_BasoAmygall_female = geneStatsall
save(geneStats_BasoAmygall_female, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_females_threeGroup.rda")

write.csv(geneStats_BasoAmygall_female, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_females_threeGroup.csv.gz"))












## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()