#####Stratify by drugs of abuse analysis#######

#drugs of abuse




#cocaine
#benzoylecgonine
#ethanol
#amphetamines
#pcp
#benzodiazepines
#barbiturates
#cocaethylene
#methanol_methyl_etoh
#isopropanol_isopropyl_etoh
#acetone_acetyl_etoh
#other_drugs
#opioids_tox
#nicotine
#cotinine
#X11_hydroxy_delta_9_thc_active
#delta_9_thc_active
#delta_9_carboxy_thc_inactive
#hallucinogens
#ketamine

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')
demo <- read.csv('csvs/rse_gene_PTSD_cohort_Feb2019.csv',header=TRUE)

#make all NA = "FALSE"
demo$benzodiazepines[is.na(demo$benzodiazepines)] <- FALSE
demo$barbiturates[is.na(demo$barbiturates)] <- FALSE
demo$cocaethylene[is.na(demo$cocaethylene)] <- FALSE
demo$acetone_acetyl_etoh[is.na(demo$acetone_acetyl_etoh)] <- FALSE
demo$other_drugs[is.na(demo$other_drugs)] <- FALSE
demo$nicotine[is.na(demo$nicotine)] <- FALSE
demo$cotinine[is.na(demo$cotinine)] <- FALSE
demo$X11_hydroxy_delta_9_thc_active[is.na(demo$X11_hydroxy_delta_9_thc_active)] <- FALSE
demo$delta_9_thc_active[is.na(demo$delta_9_thc_active)] <- FALSE
demo$delta_9_carboxy_thc_inactive[is.na(demo$delta_9_carboxy_thc_inactive)] <- FALSE
demo$hallucinogens[is.na(demo$hallucinogens)] <- FALSE
demo$ketamine[is.na(demo$ketamine)] <- FALSE


demo$any_drugs <- ifelse(demo$cocaine == "TRUE" | demo$benzoylecgonine == "TRUE" | demo$ethanol == "TRUE" | demo$amphetamines == "TRUE" | demo$pcp == "TRUE" | demo$benzodiazepines == "TRUE" | 
demo$barbiturates == "TRUE" | demo$cocaethylene == "TRUE" | demo$methanol_methyl_etoh == "TRUE" | demo$isopropanol_isopropyl_etoh == "TRUE" | demo$acetone_acetyl_etoh == "TRUE" | demo$other_drugs == "TRUE" | 
demo$opioids_tox == "1" | demo$nicotine == "TRUE" | demo$cotinine == "TRUE" | demo$X11_hydroxy_delta_9_thc_active == "TRUE" | demo$delta_9_thc_active == "TRUE" | demo$delta_9_carboxy_thc_inactive == "TRUE" | 
demo$hallucinogens == "TRUE" | demo$ketamine == "TRUE", 1, 0)


#This doesn't leave many samples with no drugs. May just try opioids_tox in future to get more balanced.


############# Drugs ##############

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

#make all NA = "FALSE"
demo$benzodiazepines[is.na(demo$benzodiazepines)] <- FALSE
demo$barbiturates[is.na(demo$barbiturates)] <- FALSE
demo$cocaethylene[is.na(demo$cocaethylene)] <- FALSE
demo$acetone_acetyl_etoh[is.na(demo$acetone_acetyl_etoh)] <- FALSE
demo$other_drugs[is.na(demo$other_drugs)] <- FALSE
demo$nicotine[is.na(demo$nicotine)] <- FALSE
demo$cotinine[is.na(demo$cotinine)] <- FALSE
demo$X11_hydroxy_delta_9_thc_active[is.na(demo$X11_hydroxy_delta_9_thc_active)] <- FALSE
demo$delta_9_thc_active[is.na(demo$delta_9_thc_active)] <- FALSE
demo$delta_9_carboxy_thc_inactive[is.na(demo$delta_9_carboxy_thc_inactive)] <- FALSE
demo$hallucinogens[is.na(demo$hallucinogens)] <- FALSE
demo$ketamine[is.na(demo$ketamine)] <- FALSE


demo$any_drugs <- ifelse(demo$cocaine == "TRUE" | demo$benzoylecgonine == "TRUE" | demo$ethanol == "TRUE" | demo$amphetamines == "TRUE" | demo$pcp == "TRUE" | demo$benzodiazepines == "TRUE" | 
demo$barbiturates == "TRUE" | demo$cocaethylene == "TRUE" | demo$methanol_methyl_etoh == "TRUE" | demo$isopropanol_isopropyl_etoh == "TRUE" | demo$acetone_acetyl_etoh == "TRUE" | demo$other_drugs == "TRUE" | 
demo$opioids_tox == "1" | demo$nicotine == "TRUE" | demo$cotinine == "TRUE" | demo$X11_hydroxy_delta_9_thc_active == "TRUE" | demo$delta_9_thc_active == "TRUE" | demo$delta_9_carboxy_thc_inactive == "TRUE" | 
demo$hallucinogens == "TRUE" | demo$ketamine == "TRUE", 1, 0)

#merge any_drugs column of demo with rse_object
colData(rse_gene) = cbind(colData(rse_gene) , demo[rse_gene$BrNum,195, drop=FALSE])


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


## also add snpPCs based on association with group
mod = model.matrix(~any_drugs + Group + Sex + AgeDeath + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate + ERCCsumLogErr + snpPC1 + snpPC2 + snpPC3 + snpPC8 + snpPC9 + snpPC10, data = colData(rse_gene))

#load cov_rse object from PTSD data
load("rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

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

modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_qsvs_anydrugs.Rdata')



#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
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
colnames(sigGene) = paste0(colnames(sigGene), "_anydrugs")






###Stuff we want

###### GENE ########


#All columns
geneStatsall = cbind(sigGene, rowData(rse_gene))
## write out qSVA-based stats
geneStats_BasoAmygall_anydrugs = geneStatsall
save(geneStats_BasoAmygall_anydrugs, file = "rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_anydrugs.rda")

write.csv(geneStats_BasoAmygall_anydrugs, file = gzfile("csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_anydrugs.csv.gz"))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

q()