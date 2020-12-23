###
library(jaffelab)
library(clusterProfiler)
library(dplyr)
library(fields)
library(lattice)

## create folder
dir.create("figures")
dir.create("tables")

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

############################
### load snRNA-seq data ####
############################

##  AMYGDALA
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/zForAnJa_AMY-sn-level-markerStats_MNT.rda",verbose=TRUE)
ts.amy = ts.amy[rownames(ts.amy) %in% deg_stats$ensemblID,]
logFDRs.amy = logFDRs.amy[rownames(logFDRs.amy) %in% deg_stats$ensemblID,]

cellSets_amy = mapply(function(t,fdr) {
	rownames(ts.amy)[t > 0 & fdr < 0.01]
}, as.data.frame(ts.amy), as.data.frame(10^logFDRs.amy))
lengths(cellSets_amy)

cellSets_amy_top2000 = lapply(as.data.frame(ts.amy), function(x) 
	rownames(ts.amy)[order(x,decreasing=TRUE)[1:2000]])

## DLPFC
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/zForAnJa_DLPFC-sn-level-markerStats_MNT.rda
",verbose=TRUE)
ts.dlpfc = ts.dlpfc[rownames(ts.dlpfc) %in% deg_stats$ensemblID,]
colnames(ts.dlpfc) = gsub("dlPFC_dlPFC_","",colnames(ts.dlpfc))

logFDRs.dlpfc = logFDRs.dlpfc[rownames(logFDRs.dlpfc) %in% deg_stats$ensemblID,]
colnames(logFDRs.dlpfc) = gsub("dlPFC_","",colnames(logFDRs.dlpfc))

cellSets_dlpfc = mapply(function(t,fdr) {
	rownames(ts.dlpfc)[t > 0 & fdr < 0.01]
}, as.data.frame(ts.dlpfc), as.data.frame(10^logFDRs.dlpfc))
lengths(cellSets_dlpfc)

cellSets_dlpfc_top2000 = lapply(as.data.frame(ts.dlpfc), function(x) 
	rownames(ts.dlpfc)[order(x,decreasing=TRUE)[1:2000]])


###############################
## input gene sets for PTSD ###
###############################

geneSets_cortex = lapply(with(deg_stats, 
	list(Cortex_either = Cortex_PValue_PTSD < 0.005,
		Cortex_up = Cortex_PValue_PTSD < 0.005 & Cortex_t_PTSD > 0,
		Cortex_down = Cortex_PValue_PTSD < 0.005 & Cortex_t_PTSD < 0,
		dlPFC_either = dlPFC_PValue_PTSD < 0.005,
		dlPFC_up = dlPFC_PValue_PTSD < 0.005 & dlPFC_t_PTSD > 0,
		dlPFC_down = dlPFC_PValue_PTSD < 0.005 & dlPFC_t_PTSD < 0,	
		dACC_either = dACC_PValue_PTSD < 0.005,
		dACC_up = dACC_PValue_PTSD < 0.005 & dACC_t_PTSD > 0,
		dACC_down = dACC_PValue_PTSD < 0.005 & dACC_t_PTSD < 0)),
			function(ii) deg_stats$ensemblID[ii])
lengths(geneSets_cortex)		
			
geneSets_amygdala = lapply(with(deg_stats, 
	list(Amygdala_either = Amygdala_PValue_PTSD < 0.005,
		Amygdala_up = Amygdala_PValue_PTSD < 0.005 & Amygdala_t_PTSD > 0,
		Amygdala_down = Amygdala_PValue_PTSD < 0.005 & Amygdala_t_PTSD < 0,
		BLA_either = BLA_PValue_PTSD < 0.005,
		BLA_up = BLA_PValue_PTSD < 0.005 & BLA_t_PTSD > 0,
		BLA_down = BLA_PValue_PTSD < 0.005 & BLA_t_PTSD < 0,	
		MeA_either = MeA_PValue_PTSD < 0.005,
		MeA_up = MeA_PValue_PTSD < 0.005 & MeA_t_PTSD > 0,
		MeA_down = MeA_PValue_PTSD < 0.005 & MeA_t_PTSD < 0)),
			function(ii) deg_stats$ensemblID[ii])
lengths(geneSets_amygdala)		
			
########################
#### do enrichment #####
########################

#######################
## any significance ###
#######################

##  cortex
enrich_stat_cortex = cellSets_dlpfc
for (i in seq(along = enrich_stat_cortex)) {
    type = deg_stats$ensemblID %in% cellSets_dlpfc[[i]]
	tabList = mclapply(geneSets_cortex, function(g) {
        tt = table(Set = factor(deg_stats$ensemblID %in% g, c(FALSE, TRUE)),
            CellType = factor(type, c(FALSE, TRUE)))
    }, mc.cores = 8)
	enrichList = lapply(tabList,fisher.test)
	
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = gsub(".odds ratio", "", rownames(o))
    enrich_stat_cortex[[i]] = o
}

enrichTab_cortex = do.call("cbind", enrich_stat_cortex)

##  amygdala
enrich_stat_amygdala = cellSets_amy
for (i in seq(along = enrich_stat_amygdala)) {
    type = deg_stats$ensemblID %in% cellSets_amy[[i]]
	tabList = mclapply(geneSets_amygdala, function(g) {
        tt = table(Set = factor(deg_stats$ensemblID %in% g, c(FALSE, TRUE)),
            CellType = factor(type, c(FALSE, TRUE)))
    }, mc.cores = 8)
	enrichList = lapply(tabList,fisher.test)
	
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = gsub(".odds ratio", "", rownames(o))
    enrich_stat_amygdala[[i]] = o
}

enrichTab_amy = do.call("cbind", enrich_stat_amygdala)

#######################
## top 2000 cell type #
#######################

##  cortex
enrich_stat_cortex_top = cellSets_dlpfc_top2000
for (i in seq(along = enrich_stat_cortex_top)) {
    type = deg_stats$ensemblID %in% cellSets_dlpfc_top2000[[i]]
	tabList = mclapply(geneSets_cortex, function(g) {
        tt = table(Set = factor(deg_stats$ensemblID %in% g, c(FALSE, TRUE)),
            CellType = factor(type, c(FALSE, TRUE)))
    }, mc.cores = 8)
	enrichList = lapply(tabList,fisher.test)
	
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = gsub(".odds ratio", "", rownames(o))
    enrich_stat_cortex_top[[i]] = o
}

enrichTab_cortex_top = do.call("cbind", enrich_stat_cortex_top)

##  amygdala
enrich_stat_amygdala_top = cellSets_amy_top2000
for (i in seq(along = enrich_stat_amygdala_top)) {
    type = deg_stats$ensemblID %in% cellSets_amy_top2000[[i]]
	tabList = mclapply(geneSets_amygdala, function(g) {
        tt = table(Set = factor(deg_stats$ensemblID %in% g, c(FALSE, TRUE)),
            CellType = factor(type, c(FALSE, TRUE)))
    }, mc.cores = 8)
	enrichList = lapply(tabList,fisher.test)
	
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = gsub(".odds ratio", "", rownames(o))
    enrich_stat_amygdala_top[[i]] = o
}

enrichTab_amy_top = do.call("cbind", enrich_stat_amygdala_top)

#########################
## look at enrichment ###
#########################
options(width=130)
pMat_cortex = enrichTab_cortex[, grep("Pval", colnames(enrichTab_cortex))]
orMat_cortex = enrichTab_cortex[, grep("OR", colnames(enrichTab_cortex))]
colnames(pMat_cortex) = gsub( ".Pval", "", colnames(pMat_cortex),fixed=TRUE)
colnames(orMat_cortex) = gsub( ".OR", "", colnames(orMat_cortex),fixed=TRUE)

t(round(-log10(pMat_cortex),1))
t(round(orMat_cortex,1))

pMat_cortex_top = enrichTab_cortex_top[, grep("Pval", colnames(enrichTab_cortex_top))]
orMat_cortex_top = enrichTab_cortex_top[, grep("OR", colnames(enrichTab_cortex_top))]
colnames(pMat_cortex_top) = gsub( ".Pval", "", colnames(pMat_cortex_top),fixed=TRUE)
colnames(orMat_cortex_top) = gsub( ".OR", "", colnames(orMat_cortex_top),fixed=TRUE)

t(round(-log10(pMat_cortex_top),1))
t(round(orMat_cortex_top,1))

###### amygdala
pMat_amy = enrichTab_amy[, grep("Pval", colnames(enrichTab_amy))]
orMat_amy = enrichTab_amy[, grep("OR", colnames(enrichTab_amy))]
colnames(pMat_amy) = gsub( ".Pval", "", colnames(pMat_amy),fixed=TRUE)
colnames(orMat_amy) = gsub( ".OR", "", colnames(orMat_amy),fixed=TRUE)

t(round(-log10(pMat_amy),1))
t(round(orMat_amy,1))

pMat_amy_top = enrichTab_amy_top[, grep("Pval", colnames(enrichTab_amy_top))]
orMat_amy_top = enrichTab_amy_top[, grep("OR", colnames(enrichTab_amy_top))]
colnames(pMat_amy_top) = gsub( ".Pval", "", colnames(pMat_amy_top),fixed=TRUE)
colnames(orMat_amy_top) = gsub( ".OR", "", colnames(orMat_amy_top),fixed=TRUE)
pMat_amy_top
t(round(-log10(pMat_amy_top),1))
t(round(orMat_amy_top,1))
