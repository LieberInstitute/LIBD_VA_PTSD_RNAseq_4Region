#Load up

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load('rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])

rse_exon$CatRegion <- ifelse(rse_exon$Region == "BasoAmyg" | rse_exon$Region == "MedialAmyg", "Amyg","Cortex")

#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]

with(colData(rse_exon), table(CatRegion))


ekeepIndex = which(rse_exon$CatRegion == "Cortex")
rse_exon <- rse_exon[, ekeepIndex]

with(colData(rse_exon), table(CatRegion))

table(rse_exon$CatRegion,rse_exon$Region)

load('rdas/PTSD_qsvs_Regionintxn_Cortex.Rdata',verbose=TRUE)

##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
        genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)

#load duplicate correlation
load("rdas/exonLevel_duplicateCorrelation_limma_forDE_PTSD.rda",verbose=TRUE)

# and then fit
fitExon = lmFit(vExon, modQsva,
        correlation=exon_dupCorr$consensus.correlation,
        block=colData(rse_exon)$BrNum)
eBExon = eBayes(fitExon)

#MDD vs controls (analysis 1)
sigExonMDD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonMDD) = paste0(colnames(sigExonMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigExonPTSD = topTable(eBExon,coef=3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonPTSD) = paste0(colnames(sigExonPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigExonDx = topTable(eBExon,coef=2:3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonDx) = paste0(colnames(sigExonDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostExon = topTable(eBayes(contrasts.fit(fitExon, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_exon))
colnames(PTSDvsMDDPostExon) = paste0(colnames(PTSDvsMDDPostExon), "_PTSDvsMDD")

#Region interaction MDD
outExon_interactionEffect_MDD = topTable(eBExon,coef=19,
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_MDD) = paste0(colnames(outExon_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outExon_interactionEffect_PTSD = topTable(eBExon,coef=20,
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_PTSD) = paste0(colnames(outExon_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigExonRegion = topTable(eBExon,coef=4,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonRegion) = paste0(colnames(sigExonRegion), "_Region")

##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
        genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)

#load duplicate correlation
load("rdas/exonLevel_duplicateCorrelation_limma_forDE_PTSD.rda",verbose=TRUE)

# and then fit
fitExon = lmFit(vExon, modQsva,
        correlation=exon_dupCorr$consensus.correlation,
        block=colData(rse_exon)$BrNum)
eBExon = eBayes(fitExon)

#MDD vs controls (analysis 1)
sigExonMDD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonMDD) = paste0(colnames(sigExonMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigExonPTSD = topTable(eBExon,coef=3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonPTSD) = paste0(colnames(sigExonPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigExonDx = topTable(eBExon,coef=2:3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonDx) = paste0(colnames(sigExonDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostExon = topTable(eBayes(contrasts.fit(fitExon, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_exon))
colnames(PTSDvsMDDPostExon) = paste0(colnames(PTSDvsMDDPostExon), "_PTSDvsMDD")

#Region interaction MDD
outExon_interactionEffect_MDD = topTable(eBExon,coef=19,
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_MDD) = paste0(colnames(outExon_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outExon_interactionEffect_PTSD = topTable(eBExon,coef=20,
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_PTSD) = paste0(colnames(outExon_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigExonRegion = topTable(eBExon,coef=4,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonRegion) = paste0(colnames(sigExonRegion), "_Region")


#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

load('rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])

rse_exon$onlyPTSD = ifelse(rse_exon$Group == "PTSD" ,1, 0) 

rse_exon$CatRegion <- ifelse(rse_exon$Region == "BasoAmyg" | rse_exon$Region == "MedialAmyg", "Amyg","Cortex")

#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]

ekeepIndex = which(rse_exon$CatRegion == "Cortex")
rse_exon <- rse_exon[, ekeepIndex]

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD_Cortex.Rdata',verbose=TRUE)
##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
        genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)

#load duplicate correlation
load("rdas/exonLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda",verbose=TRUE)

# and then fit
fitExon = lmFit(vExon, modQsva,
        correlation=exon_dupCorr$consensus.correlation,
        block=colData(rse_exon)$BrNum)
eBExon = eBayes(fitExon)

#PTSD only
sigExononlyPTSD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExononlyPTSD) = paste0(colnames(sigExononlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outExon_interactionEffect_onlyPTSD = topTable(eBExon,coef=18,
        p.value = 1,number=nrow(rse_exon), sort="none")
colnames(outExon_interactionEffect_onlyPTSD) = paste0(colnames(outExon_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")

#Region effect
sigExonRegion_onlyPTSD = topTable(eBExon,coef=3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonRegion_onlyPTSD) = paste0(colnames(sigExonRegion_onlyPTSD), "_Region_onlyPTSD")

###### EXON ########
#All columns
exonStatsall = cbind(sigExonPTSD, sigExonMDD, sigExonDx, sigExononlyPTSD, PTSDvsMDDPostExon, sigExonRegion, outExon_interactionEffect_MDD, 
	outExon_interactionEffect_PTSD, outExon_interactionEffect_onlyPTSD,sigExonRegion_onlyPTSD)
exonStatsall = cbind(exonStatsall, rowData(rse_exon))
save(exonStatsall, file = "rdas/all_regions/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda")

library(utils)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()

