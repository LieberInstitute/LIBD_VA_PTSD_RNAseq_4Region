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

load('rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])

rse_tx$CatRegion <- ifelse(rse_tx$Region == "BasoAmyg" | rse_tx$Region == "MedialAmyg", "Amyg","Cortex")

#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]

with(colData(rse_tx), table(CatRegion))

tkeepIndex = which(rse_tx$CatRegion == "Cortex")
rse_tx <- rse_tx[, tkeepIndex]

with(colData(rse_tx), table(CatRegion))

table(rse_tx$CatRegion,rse_tx$Region)

load('rdas/PTSD_qsvs_Regionintxn_Cortex.Rdata',verbose=TRUE)

##### TRANSCRIPT ######

txExprs = log2(assays(rse_tx)$tpm+ 1)

#load duplicate correlation
load("rdas/txLevel_duplicateCorrelation_limma_forDE_PTSD.rda",verbose=TRUE)

# and then fit
fitTx = lmFit(txExprs, modQsva,
        correlation=tx_dupCorr$consensus.correlation,
        block=colData(rse_tx)$BrNum)
eBTx = eBayes(fitTx)

#MDD vs controls (analysis 1)
sigTxMDD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxMDD) = paste0(colnames(sigTxMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigTxPTSD = topTable(eBTx,coef=3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxPTSD) = paste0(colnames(sigTxPTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigTxDx = topTable(eBTx,coef=2:3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxDx) = paste0(colnames(sigTxDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostTx = topTable(eBayes(contrasts.fit(fitTx, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, genelist = rowRanges(rse_tx), sort="none", n = nrow(rse_tx))
colnames(PTSDvsMDDPostTx) = paste0(colnames(PTSDvsMDDPostTx), "_PTSDvsMDD")

#Region interaction MDD
outTx_interactionEffect_MDD = topTable(eBTx,coef=19,
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_MDD) = paste0(colnames(outTx_interactionEffect_MDD), "_rintxnMDD")
												 
#Region interaction PTSD
outTx_interactionEffect_PTSD = topTable(eBTx,coef=20,
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_PTSD) = paste0(colnames(outTx_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigTxRegion = topTable(eBTx,coef=4,
	p.value = 1,number=nrow(rse_tx), sort="none")
colnames(sigTxRegion) = paste0(colnames(sigTxRegion), "_Region")

#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

load('rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])

rse_tx$onlyPTSD = ifelse(rse_tx$Group == "PTSD" ,1, 0) 

rse_tx$CatRegion <- ifelse(rse_tx$Region == "BasoAmyg" | rse_tx$Region == "MedialAmyg", "Amyg","Cortex")

#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]

tkeepIndex = which(rse_tx$CatRegion == "Cortex")
rse_tx <- rse_tx[, tkeepIndex]

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD_Cortex.Rdata',verbose=TRUE)

##### TRANSCRIPT ######

txExprs = log2(assays(rse_tx)$tpm+ 1)

#load duplicate correlation
load("rdas/txLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda",verbose=TRUE)

# and then fit
fitTx = lmFit(txExprs, modQsva,
        correlation=tx_dupCorr$consensus.correlation,
        block=colData(rse_tx)$BrNum)
eBTx = eBayes(fitTx)

#PTSD only
sigTxonlyPTSD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxonlyPTSD) = paste0(colnames(sigTxonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outTx_interactionEffect_onlyPTSD = topTable(eBTx,coef=18,
        p.value = 1,number=nrow(rse_tx), sort="none")
colnames(outTx_interactionEffect_onlyPTSD) = paste0(colnames(outTx_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")
								
#Region effect
sigTxRegion_onlyPTSD = topTable(eBTx,coef=3,
	p.value = 1,number=nrow(rse_tx), sort="none")
colnames(sigTxRegion_onlyPTSD) = paste0(colnames(sigTxRegion_onlyPTSD), "_Region_onlyPTSD")

####### TRANSCRIPT ########
#All columns
txStatsall = cbind(sigTxPTSD, sigTxMDD, sigTxDx, sigTxonlyPTSD, PTSDvsMDDPostTx, sigTxRegion, outTx_interactionEffect_MDD, 
	outTx_interactionEffect_PTSD, outTx_interactionEffect_onlyPTSD,sigTxRegion_onlyPTSD)
txStatsall = cbind(txStatsall, rowData(rse_tx))
save(txStatsall, file = "rdas/all_regions/txStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda")

library(utils)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()

