#####################################
#Make qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC.R
#and https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC_allFeatures.R
#and /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/de_analysis.R)
#See also '.../R_scripts/qSV_model_analysis_dACC_PTSD_updated.R completed first
#####################################

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
load('rdas/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
load('rdas/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])
#exon
colData(rse_exon) = cbind(colData(rse_exon) , mds[rse_exon$BrNum,])
#junction
colData(rse_jx) = cbind(colData(rse_jx) , mds[rse_jx$BrNum,])
#transcript
colData(rse_tx) = cbind(colData(rse_tx) , mds[rse_tx$BrNum,])


## expression filter to remove lowly expressed stuff 
## do across all regions so we're looking at the same features
####Do not filter for 23andMe -- need to filter later for other analyses though

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]
#exon
eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.2
rse_exon <- rse_exon[eIndex , ]
#junction
rowRanges(rse_jx)$Length <- 100
jIndex = rowMeans(getRPKM(rse_jx, "Length")) > 0.75 & rowData(rse_jx)$Class != "Novel"
rse_jx <- rse_jx[jIndex , ]
#transcript
tIndex = rowMeans(assays(rse_tx)$tpm) > 0.2
rse_tx <- rse_tx[tIndex , ]


#Load QSVA
load('rdas/PTSD_qsvs.Rdata')


#filter for dACC
#gene
keepIndex = which(rse_gene$Region == "dACC")
rse_gene <- rse_gene[, keepIndex]
#exon
rse_exon <- rse_exon[, keepIndex]
#junction
rse_jx <- rse_jx[, keepIndex]
#transcript
rse_tx <- rse_tx[, keepIndex]

#Get rid of region columns
colIndex <-  !grepl("Region", colnames(mod))
mod <- mod[keepIndex, colIndex]
modQsva <- modQsva[keepIndex,!grepl("Region", colnames(modQsva))]
colnames(modQsva)[c(1, 18:36)] = c("Int", gsub("PC", "qSV", colnames(modQsva)[18:36]))



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



##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)
fitExon = lmFit(vExon)
eBExon = eBayes(fitExon)

#MDD vs controls (analysis 1)
sigExonMDD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonMDD) = paste0(colnames(sigExonMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigExonPTSD = topTable(eBExon,coef=3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonPTSD) = paste0(colnames(sigExonPTSD), "_PTSD")

#Interaction (analysis 3)
sigExonDx = topTable(eBExon,coef=2:3,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExonDx) = paste0(colnames(sigExonDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostExon = topTable(eBayes(contrasts.fit(fitExon, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_exon))
colnames(PTSDvsMDDPostExon) = paste0(colnames(PTSDvsMDDPostExon), "_PTSDvsMDD")


##### JUNCTION ######



dje = DGEList(counts = assays(rse_jx)$counts,
	genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)
fitJxn = lmFit(vJxn)
eBJxn = eBayes(fitJxn)


#MDD vs controls (analysis 1)
sigJxnMDD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnMDD) = paste0(colnames(sigJxnMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigJxnPTSD = topTable(eBJxn,coef=3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnPTSD) = paste0(colnames(sigJxnPTSD), "_PTSD")

#Interaction (analysis 3)
sigJxnDx = topTable(eBJxn,coef=2:3,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnDx) = paste0(colnames(sigJxnDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostJxn = topTable(eBayes(contrasts.fit(fitJxn, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, sort="none", n = nrow(rse_jx))
colnames(PTSDvsMDDPostJxn) = paste0(colnames(PTSDvsMDDPostJxn), "_PTSDvsMDD")



##### TRANSCRIPT ######

fitTx = lmFit(log2(assays(rse_tx)$tpm + 1), modQsva)
eBTx = eBayes(fitTx)


#MDD vs controls (analysis 1)
sigTxMDD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxMDD) = paste0(colnames(sigTxMDD), "_MDD")

#PTSD vs controls (analysis 1)
sigTxPTSD = topTable(eBTx,coef=3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxPTSD) = paste0(colnames(sigTxPTSD), "_PTSD")

#Interaction (analysis 3)
sigTxDx = topTable(eBTx,coef=2:3,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxDx) = paste0(colnames(sigTxDx), "_ANOVA")

#PTSDvsMDD using limma
PTSDvsMDDContrast <- makeContrasts(GroupPTSD-GroupMDD,levels=modQsva)
PTSDvsMDDPostTx = topTable(eBayes(contrasts.fit(fitTx, PTSDvsMDDContrast)),
    coef=1,  p.value = 1, genelist = rowRanges(rse_tx), sort="none", n = nrow(rse_tx))
colnames(PTSDvsMDDPostTx) = paste0(colnames(PTSDvsMDDPostTx), "_PTSDvsMDD")



#onlyPTSD analysis (analysis 2)
# 1. redefine model
# 2. add to geneStats


#Define by PTSD 
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 
rse_exon$onlyPTSD = ifelse(rse_exon$Group == "PTSD" ,1, 0) 
rse_jx$onlyPTSD = ifelse(rse_jx$Group == "PTSD" ,1, 0) 
rse_tx$onlyPTSD = ifelse(rse_tx$Group == "PTSD" ,1, 0) 

load('rdas/PTSD_onlyPTSD_qsvs.Rdata')


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


##### EXON ######

dee = DGEList(counts = assays(rse_exon)$counts,
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)
fitExon = lmFit(vExon)
eBExon = eBayes(fitExon)

#PTSD only
sigExononlyPTSD = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon), sort="none")
colnames(sigExononlyPTSD) = paste0(colnames(sigExononlyPTSD), "_onlyPTSD")


##### JUNCTION ######

dje = DGEList(counts = assays(rse_jx)$counts,
	genes = rowData(rse_jx))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)
fitJxn = lmFit(vJxn)
eBJxn = eBayes(fitJxn)


#PTSD only
sigJxnonlyPTSD = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jx), sort="none")
colnames(sigJxnonlyPTSD) = paste0(colnames(sigJxnonlyPTSD), "_onlyPTSD")


##### TRANSCRIPT ######

fitTx = lmFit(log2(assays(rse_tx)$tpm + 1), modQsva)
eBTx = eBayes(fitTx)


#PTSD only
sigTxonlyPTSD = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx), genelist = rowRanges(rse_tx), sort="none")
colnames(sigTxonlyPTSD) = paste0(colnames(sigTxonlyPTSD), "_onlyPTSD")



###Stuff we want

###### GENE ########

#Merge qSVA into one big genestats dataframe
geneStats = cbind(sigGenePTSD[,c(11,13:15)], sigGeneMDD[,c(11,13:15)],sigGeneDx[,c(14:16)],sigGeneonlyPTSD[,c(11,13:15)], PTSDvsMDDPostGene[,c(11,13:15)])
geneStats = cbind(geneStats, rowData(rse_gene))


## write out qSVA-based stats
geneStats_dACC = geneStats
save(geneStats_dACC, file = "rdas/dACC/geneStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_dACC, file = gzfile("csvs/dACC/geneStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))


#All columns
geneStatsall = cbind(sigGenePTSD, sigGeneMDD, sigGeneDx, sigGeneonlyPTSD, PTSDvsMDDPostGene)
geneStatsall = cbind(geneStatsall, rowData(rse_gene))
## write out qSVA-based stats
geneStats_dACCall = geneStatsall
save(geneStats_dACCall, file = "rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(geneStats_dACCall, file = gzfile("csvs/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))



###### EXON ########

#Merge qSVA into one big genestats dataframe
exonStats = cbind(sigExonPTSD[,c(16,18:20)], sigExonMDD[,c(16,18:20)],sigExonDx[,c(19:21)],sigExononlyPTSD[,c(16,18:20)], PTSDvsMDDPostExon[,c(16,18:20)])
exonStats = cbind(exonStats, rowData(rse_exon))


## write out qSVA-based stats
exonStats_dACC = exonStats
save(exonStats_dACC, file = "rdas/dACC/exonStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(exonStats_dACC, file = gzfile("csvs/dACC/exonStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))


#All columns
exonStatsall = cbind(sigExonPTSD, sigExonMDD, sigExonDx, sigExononlyPTSD, PTSDvsMDDPostExon)
exonStatsall = cbind(exonStatsall, rowData(rse_exon))
## write out qSVA-based stats
exonStats_dACCall = exonStatsall
save(exonStats_dACCall, file = "rdas/dACC/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(exonStats_dACCall, file = gzfile("csvs/dACC/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))



###### JUNCTION ########

#Merge qSVA into one big genestats dataframe
jxnStats = cbind(sigJxnPTSD[,c(17,19:21)], sigJxnMDD[,c(17,19:21)],sigJxnDx[,c(20:22)],sigJxnonlyPTSD[,c(17,19:21)], PTSDvsMDDPostJxn[,c(17,19:21)])
jxnStats = cbind(jxnStats, rowData(rse_jx))


## write out qSVA-based stats
jxnStats_dACC = jxnStats
save(jxnStats_dACC, file = "rdas/dACC/jxStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(jxnStats_dACC, file = gzfile("csvs/dACC/jxStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))

#All columns
jxnStatsall = cbind(sigJxnPTSD, sigJxnMDD, sigJxnDx, sigJxnonlyPTSD, PTSDvsMDDPostJxn)
jxnStatsall = cbind(jxnStatsall, rowData(rse_jx))
## write out qSVA-based stats
jxnStats_dACCall = jxnStatsall
save(jxnStats_dACCall, file = "rdas/dACC/jxStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(jxnStats_dACCall, file = gzfile("csvs/dACC/jxStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))

####### TRANSCRIPT ########

#Merge qSVA into one big genestats dataframe
txStats = cbind(sigTxPTSD[,c(28,30:32)], sigTxMDD[,c(28,30:32)],sigTxDx[,c(31:33)],sigTxonlyPTSD[,c(28,30:32)], PTSDvsMDDPostTx[,c(28,30:32)])
txStats = cbind(txStats, rowData(rse_tx))


## write out qSVA-based stats
txStats_dACC = txStats
save(txStats_dACC, file = "rdas/dACC/txStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(txStats_dACC, file = gzfile("csvs/dACC/txStats_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))


#All columns
txStatsall = cbind(sigTxPTSD, sigTxMDD, sigTxDx, sigTxonlyPTSD, PTSDvsMDDPostTx)
txStatsall = cbind(txStatsall, rowData(rse_tx))
## write out qSVA-based stats
txStats_dACCall = txStatsall
save(txStats_dACCall, file = "rdas/dACC/txStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda")
#dir.create("csvs",showWarnings = FALSE)
write.csv(txStats_dACCall, file = gzfile("csvs/dACC/txStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.csv.gz"))


q()