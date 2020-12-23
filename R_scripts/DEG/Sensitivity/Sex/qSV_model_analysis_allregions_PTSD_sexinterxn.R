#load up

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)


setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

###BasoAmyg

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])


## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]


load('rdas/PTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)


###Analysis

#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


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

#Diagnosis interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")

#Sex main effect
sigGeneSex = topTable(eBGene,coef=4,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneSex) = paste0(colnames(sigGeneSex), "_Sex")

#Sex interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=18,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_sexintxnMDD")

#Sex interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=19,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_sexintxnPTSD")




#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse object from PTSD data
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata', verbose = TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)

#filter for BasoAmyg 
keepIndex = which(rse_gene$Region == "BasoAmyg")
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
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")

#Sex interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=17,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_sexintxnonlyPTSD")




###Stuff we want

###### GENE ########

#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPost, sigGeneSex, outGene_interactionEffect_MDD, outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.rda")
write.csv(geneStatsall, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.csv.gz"))

##########################################################

#MedialAmyg

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

load('rdas/PTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)


###Analysis

#filter for MedialAmyg 
keepIndex = which(rse_gene$Region == "MedialAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


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

#Diagnosis interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")

#Sex main effect
sigGeneSex = topTable(eBGene,coef=4,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneSex) = paste0(colnames(sigGeneSex), "_Sex")


#Sex interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=18,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_sexintxnMDD")

#Sex interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=19,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_sexintxnPTSD")



#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)

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
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


#Sex interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=17,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_sexintxnonlyPTSD")


###Stuff we want

###### GENE ########

#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPost, sigGeneSex, outGene_interactionEffect_MDD, outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/MedialAmyg/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.rda")
write.csv(geneStatsall, file = gzfile("csvs/MedialAmyg/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.csv.gz"))

##########################################################

#dACC

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

load('rdas/PTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)


###Analysis

#filter for dACC 
keepIndex = which(rse_gene$Region == "dACC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


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

#Diagnosis interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")

#Sex main effect
sigGeneSex = topTable(eBGene,coef=4,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneSex) = paste0(colnames(sigGeneSex), "_Sex")

#Sex interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=18,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_sexintxnMDD")

#Sex interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=19,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_sexintxnPTSD")



#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)

#filter for dACC 
keepIndex = which(rse_gene$Region == "dACC")
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
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


#Sex interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=17,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_sexintxnonlyPTSD")


###Stuff we want

###### GENE ########

#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPost, sigGeneSex, outGene_interactionEffect_MDD, outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/dACC/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.rda")
write.csv(geneStatsall, file = gzfile("csvs/dACC/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.csv.gz"))


##########################################################

#DLPFC

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

load('rdas/PTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)


###Analysis

#filter for DLPFC 
keepIndex = which(rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]


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

#Diagnosis interaction (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPost = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPost) = paste0(colnames(PTSDvsMDDPost), "_PTSDvsMDD")

#Sex main effect
sigGeneSex = topTable(eBGene,coef=4,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneSex) = paste0(colnames(sigGeneSex), "_Sex")

#Sex interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=18,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_sexintxnMDD")

#Sex interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=19,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_sexintxnPTSD")



#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#load rse_gene object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda",verbose=TRUE)
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

## expression filter to remove lowly expressed stuff (can do later as well)
##		do across all regions so we're looking at the same genes
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs_sexinterxn.Rdata',verbose=TRUE)

#filter for DLPFC 
keepIndex = which(rse_gene$Region == "DLPFC")
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
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")


#Sex interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=17,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_sexintxnonlyPTSD")


###Stuff we want

###### GENE ########

#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPost, sigGeneSex, outGene_interactionEffect_MDD, outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/DLPFC/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.rda")
write.csv(geneStatsall, file = gzfile("csvs/DLPFC/geneStats_allcols_DE_qSVA_sexinterxn_threeGroup.csv.gz"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()
