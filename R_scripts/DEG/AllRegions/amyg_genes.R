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

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

with(colData(rse_gene), table(CatRegion))

#filter for Amyg
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

with(colData(rse_gene), table(CatRegion))

table(rse_gene$CatRegion,rse_gene$Region)

load('rdas/PTSD_qsvs_Regionintxn_Amyg.Rdata',verbose=TRUE)

##### GENE #####

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=FALSE)

#load duplicate correlation
load("rdas/geneLevel_duplicateCorrelation_limma_forDE_PTSD.rda",verbose=TRUE)

# and then fit
fitGene = lmFit(vGene, modQsva,
        correlation=gene_dupCorr$consensus.correlation,
        block=colData(rse_gene)$BrNum)
eBGene = eBayes(fitGene)

#MDD vs controls (analysis 1)
sigGeneMDD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneMDD) = paste0(colnames(sigGeneMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigGenePTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGenePTSD) = paste0(colnames(sigGenePTSD), "_PTSD")

#Interaction of diagnosis (analysis 3)
sigGeneDx = topTable(eBGene,coef=2:3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneDx) = paste0(colnames(sigGeneDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostGene = topTable(eBayes(contrasts.fit(fitGene, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_gene))
colnames(PTSDvsMDDPostGene) = paste0(colnames(PTSDvsMDDPostGene), "_PTSDvsMDD")

#Region interaction MDD
outGene_interactionEffect_MDD = topTable(eBGene,coef=19,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_MDD) = paste0(colnames(outGene_interactionEffect_MDD), "_rintxnMDD")

#Region interaction PTSD
outGene_interactionEffect_PTSD = topTable(eBGene,coef=20,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_PTSD) = paste0(colnames(outGene_interactionEffect_PTSD), "_rintxnPTSD")

#Region effect
sigGeneRegion = topTable(eBGene,coef=4,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneRegion) = paste0(colnames(sigGeneRegion), "_Region")

#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats

#reload rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

#Define by PTSD
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")

## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#filter for Amyg
gkeepIndex = which(rse_gene$CatRegion == "Amyg")
rse_gene <- rse_gene[, gkeepIndex]

load('rdas/PTSD_qsvs_Regionintxn_onlyPTSD_Amyg.Rdata',verbose=TRUE)

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)

load("rdas/geneLevel_duplicateCorrelation_limma_forDE_PTSD_onlyPTSD.rda",verbose=TRUE)

# and then fit
fitGene = lmFit(vGene, modQsva,
        correlation=gene_dupCorr$consensus.correlation,
        block=colData(rse_gene)$BrNum)
eBGene = eBayes(fitGene)

#PTSD only
sigGeneonlyPTSD = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneonlyPTSD) = paste0(colnames(sigGeneonlyPTSD), "_onlyPTSD")

#Region interaction onlyPTSD
outGene_interactionEffect_onlyPTSD = topTable(eBGene,coef=18,
        p.value = 1,number=nrow(rse_gene), sort="none")
colnames(outGene_interactionEffect_onlyPTSD) = paste0(colnames(outGene_interactionEffect_onlyPTSD), "_rintxnonlyPTSD")

#Region effect
sigGeneRegion_onlyPTSD = topTable(eBGene,coef=3,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGeneRegion_onlyPTSD) = paste0(colnames(sigGeneRegion_onlyPTSD), "_Region_onlyPTSD")

###### GENE ########
#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPostGene, sigGeneRegion, outGene_interactionEffect_MDD, 
	outGene_interactionEffect_PTSD, outGene_interactionEffect_onlyPTSD,sigGeneRegion_onlyPTSD)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
save(geneStatsall, file = "rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()

