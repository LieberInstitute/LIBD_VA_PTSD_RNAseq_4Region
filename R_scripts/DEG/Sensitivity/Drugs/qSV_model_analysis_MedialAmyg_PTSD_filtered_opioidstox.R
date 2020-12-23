#####Stratify by drugs of abuse analysis#######

###opioids_tox


library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse object
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

#load demographics data
demo <- read.csv('csvs/rse_gene_PTSD_cohort_Feb2019.csv',header=TRUE)
demo$BrNum <- demo$brnum
rownames(demo) <- demo$BrNum


#merge any_drugs column of demo with rse_object
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,84, drop=FALSE])


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

#load mod, qSVs
## also add snpPCs based on association with group
load('rdas/PTSD_qsvs_opioidstox.Rdata')

#filter for MedialAmyg 
#gene
keepIndex = which(rse_gene$Region == "MedialAmyg")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[1] = "Int"



##### GENE ######

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))

#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#drugs vs no drugs
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
colnames(sigGene) = paste0(colnames(sigGene), "_opioidstox")


###Stuff we want

###### GENE ########


#All columns
geneStatsall = cbind(sigGene, rowData(rse_gene))
## write out qSVA-based stats
geneStats_MedialAmygall_opioidstox = geneStatsall
save(geneStats_MedialAmygall_opioidstox, file = "rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_opioidstox.rda")

write.csv(geneStats_MedialAmygall_opioidstox, file = gzfile("csvs/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_opioidstox.csv.gz"))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()