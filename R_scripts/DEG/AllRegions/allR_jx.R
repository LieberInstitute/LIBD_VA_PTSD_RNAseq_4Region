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


load('rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])

#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]

load('rdas/PTSD_qsvs_Regionintxn.Rdata',verbose=TRUE)

##### JUNCTION ######

dje = DGEList(counts = assays(rse_jx)$counts,
        genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)

## do duplicate correlation
#jx_dupCorr = duplicateCorrelation(vJxn$E, modQsva, block=colData(rse_jx)$BrNum)
#save(jx_dupCorr, file = "rdas/jxLevel_duplicateCorrelation_limma_forDE_PTSD.rda")
load("rdas/jxLevel_duplicateCorrelation_limma_forDE_PTSD.rda",verbose=TRUE)

# and then fit
fitJxn = lmFit(vJxn, modQsva,
        correlation=jx_dupCorr$consensus.correlation,
        block=colData(rse_jx)$BrNum)
eBJxn = eBayes(fitJxn)

#MDD vs controls (analysis 1)
sigJxnMDD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnMDD) = paste0(colnames(sigJxnMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigJxnPTSD = topTable(eBJxn,coef=3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnPTSD) = paste0(colnames(sigJxnPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigJxnDx = topTable(eBJxn,coef=2:3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnDx) = paste0(colnames(sigJxnDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostJxn = topTable(eBayes(contrasts.fit(fitJxn, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_jx))
colnames(PTSDvsMDDPostJxn) = paste0(colnames(PTSDvsMDDPostJxn), "_PTSDvsMDD")

#Region interaction MDD
outJxn_interactionEffect_MDD = topTable(eBJxn,coef=c(21,23,25),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_MDD) = paste0(colnames(outJxn_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outJxn_interactionEffect_PTSD = topTable(eBJxn,coef=c(22,24,26),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_PTSD) = paste0(colnames(outJxn_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigJxnRegion = topTable(eBJxn,coef=c(4,5,6),
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnRegion) = paste0(colnames(sigJxnRegion), "_Region")

#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

load('rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])

#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]

rse_jx$onlyPTSD = ifelse(rse_jx$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD.Rdata',verbose=TRUE)

##### JUNCTION ######

dje = DGEList(counts = assays(rse_jx)$counts,
        genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)

## do duplicate correlation
jx_dupCorr = duplicateCorrelation(vJxn$E, modQsva, block=colData(rse_jx)$BrNum)
save(jx_dupCorr, file = "rdas/jxLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda")

# and then fit
fitJxn = lmFit(vJxn, modQsva,
        correlation=jx_dupCorr$consensus.correlation,
        block=colData(rse_jx)$BrNum)
eBJxn = eBayes(fitJxn)

#PTSD only
sigJxnonlyPTSD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnonlyPTSD) = paste0(colnames(sigJxnonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outJxn_interactionEffect_onlyPTSD = topTable(eBJxn,coef=c(20:22),
        p.value = 1,number=nrow(rse_jx), sort="none")
colnames(outJxn_interactionEffect_onlyPTSD) = paste0(colnames(outJxn_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")
													  
#Region effect
sigJxnRegion_onlyPTSD = topTable(eBJxn,coef=c(3,4,5),
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnRegion_onlyPTSD) = paste0(colnames(sigJxnRegion_onlyPTSD), "_Region_onlyPTSD")

###### JUNCTION ########
#All columns
jxnStatsall = cbind(sigJxnPTSD, sigJxnMDD, sigJxnDx, sigJxnonlyPTSD, PTSDvsMDDPostJxn, sigJxnRegion, outJxn_interactionEffect_MDD, 
	outJxn_interactionEffect_PTSD, outJxn_interactionEffect_onlyPTSD,sigJxnRegion_onlyPTSD)
jxnStatsall = cbind(jxnStatsall, rowData(rse_jx))
## write out qSVA-based stats
save(jxnStatsall, file = "rdas/all_regions/jxStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda")



library(utils)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()

