######################
## WGCNA synthesis ###
######################
library(jaffelab)

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

#########################
## load wgcna stats #####
#########################
wgcna_stats_files = list(
		Cortex = "../../Data/rdas/WGCNA/AllRegions/MEvsDx_Cortex.rda",
		dlPFC = "../../Data/rdas/WGCNA/DLPFC/MEvsDx_DLPFC.rda",
		dACC = "../../Data/rdas/WGCNA/dACC/MEvsDx_dACC.rda",
		Amygdala = "../../Data/rdas/WGCNA/AllRegions/MEvsDx_Amyg.rda",
		MeA = "../../Data/rdas/WGCNA/MedialAmyg/MEvsDx_MedialAmyg.rda",
		BLA = "../../Data/rdas/WGCNA/BasoAmyg/MEvsDx_BasoAmyg.rda")
wgcna_stats_ptsd = lapply(wgcna_stats_files, function(x) {
	load(x)
	coefAdj_PTSD
})
wgcna_stats_mdd = lapply(wgcna_stats_files, function(x) {
	load(x)
	coefAdj_MDD
})
wgcna_stats_onlyptsd = lapply(wgcna_stats_files, function(x) {
	load(x)
	coefAdj_onlyPTSD
})


#### invert to regions ###
regions = names(wgcna_stats_files)
names(regions) = regions
regionList = lapply(regions, function(r) {
	ptsd = wgcna_stats_ptsd[[r]] 
	mdd = wgcna_stats_mdd[[r]] 
	onlyptsd = wgcna_stats_onlyptsd[[r]] 
	m = paste0("ME", 0:(nrow(ptsd)-1))
	tMat = data.frame(PTSD = ptsd[m,"t value"], 
		MDD = mdd[m,"t value"], onlyPTSD = onlyptsd[m,"t value"])
	colnames(tMat) = paste0(colnames(tMat), "_t")
	pMat = data.frame(PTSD = ptsd[m,"Pr(>|t|)"],
		MDD = mdd[m,"Pr(>|t|)"], onlyPTSD = onlyptsd[m,"Pr(>|t|)"])
	colnames(pMat) = paste0(colnames(pMat), "_p")
	o = cbind(tMat, pMat)
	o[,c(1,4,2,5,3,6)]
})

##merge
region_stats =do.call("rbind", regionList)
region_stats$Module_Name = gsub("\\.", "_", rownames(region_stats))


############################
### load other stuff #######

wgcna_modules_files = list(
		Cortex = "../../Data/rdas/WGCNA/AllRegions/constructed_network_signed_bicor_Cortex.rda",
		dlPFC = "../../Data/rdas/WGCNA/DLPFC/constructed_network_signed_bicor_DLPFC.rda",
		dACC = "../../Data/rdas/WGCNA/dACC/constructed_network_signed_bicor_dACC.rda",
		Amygdala = "../../Data/rdas/WGCNA/AllRegions/constructed_network_signed_bicor_Amyg.rda",
		MeA = "../../Data/rdas/WGCNA/MedialAmyg/constructed_network_signed_bicor_MedialAmyg.rda",
		BLA = "../../Data/rdas/WGCNA/BasoAmyg/constructed_network_signed_bicor_BasoAmyg.rda")
wgcna_modules_list = lapply(wgcna_modules_files, function(x) {
	load(x)
	return(net)
})

## add number of modules 
gene_numbers = unlist(lapply(wgcna_modules_list, function(x) {
	tt = table(x$colors)
	names(tt) = paste0("ME", names(tt))
	tt
}))
region_stats$numGenes = gene_numbers[rownames(region_stats)]
region_stats = region_stats[,c(7:8, 1:6)] # reorder

## get gene membership
gene_member_mat = sapply(wgcna_modules_list, function(x) x$colors)
identical(rownames(deg_stats), rownames(gene_member_mat)) # TRUE

## write out
gene_member_out = as.data.frame(gene_member_mat)
gene_member_out$Symbol = deg_stats$Symbol
gene_member_out$GencodeID = rownames(gene_member_out)
gene_member_out = gene_member_out[,c(8,7,1:6)]
write.csv(gene_member_out, "tables/suppTable_WGCNA_gene_membership.csv",row.names=FALSE)

########################### 
## enrichment of DE genes #
###########################

##### overlap among significant genes
geneSets_PTSD = with(deg_stats, 
	list(Cortex_PTSD = Cortex_PValue_PTSD < 0.005,
		dlPFC_PTSD = dlPFC_PValue_PTSD < 0.005,
		dACC_PTSD = dACC_PValue_PTSD < 0.005,
		Amygdala_PTSD = Amygdala_PValue_PTSD < 0.005,
		MeA_PTSD = MeA_PValue_PTSD < 0.005,
		BLA_PTSD = BLA_PValue_PTSD < 0.005))
geneSets_MDD = with(deg_stats, 
	list(Cortex_MDD = Cortex_PValue_MDD < 0.005,
		dlPFC_MDD = dlPFC_PValue_MDD < 0.005,
		dACC_MDD = dACC_PValue_MDD < 0.005,
		Amygdala_MDD = Amygdala_PValue_MDD < 0.005,
		MeA_MDD = MeA_PValue_MDD < 0.005,
		BLA_MDD = BLA_PValue_MDD < 0.005))
	
enrichStats_list_PTSD = mapply(function(g, de) {
	gIndexes = splitit(g) 
	tabList = lapply(gIndexes, function(ii) {
		z = rep(FALSE,nrow(deg_stats))
		z[ii] = TRUE
		tt = table(inMod = z, isDE = de)
	})
	enrichList = lapply(tabList, fisher.test)
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = paste0("ME", gsub(".odds ratio", "", rownames(o)))
	return(o)
}, as.data.frame(gene_member_mat), geneSets_PTSD,SIMPLIFY=FALSE)

enrichStats_list_MDD = mapply(function(g, de) {
	gIndexes = splitit(g) 
	tabList = lapply(gIndexes, function(ii) {
		z = rep(FALSE,nrow(deg_stats))
		z[ii] = TRUE
		tt = table(inMod = z, isDE = de)
	})
	enrichList = lapply(tabList, fisher.test)
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = paste0("ME", gsub(".odds ratio", "", rownames(o)))
	return(o)
}, as.data.frame(gene_member_mat), geneSets_MDD,SIMPLIFY=FALSE)

## merge 
enrichStats_PTSD = do.call("rbind", enrichStats_list_PTSD)
colnames(enrichStats_PTSD) = paste0("PTSD_", colnames(enrichStats_PTSD))
enrichStats_MDD = do.call("rbind", enrichStats_list_MDD)
colnames(enrichStats_MDD) = paste0("MDD_", colnames(enrichStats_MDD))
enrichStats = cbind(enrichStats_PTSD[rownames(region_stats),],
				enrichStats_MDD[rownames(region_stats),])
region_stats = cbind(region_stats, enrichStats)

#####################################
### gene set enrichment of modules ##
#####################################
library(clusterProfiler)
wgcna_go_files = list(
		Cortex = "../../Data/rdas/WGCNA/AllRegions/wgcna_Cortex_GO_clusterProfiler.rda",
		dlPFC = "../../Data/rdas/WGCNA/DLPFC/wgcna_DLPFC_GO_clusterProfiler.rda",
		dACC = "../../Data/rdas/WGCNA/dACC/wgcna_dACC_GO_clusterProfiler.rda",
		Amygdala = "../../Data/rdas/WGCNA/AllRegions/wgcna_Amyg_GO_clusterProfiler.rda",
		MeA = "../../Data/rdas/WGCNA/MedialAmyg/wgcna_MedialAmyg_GO_clusterProfiler.rda",
		BLA = "../../Data/rdas/WGCNA/BasoAmyg/wgcna_BasoAmyg_GO_clusterProfiler.rda")
## read in
wgcna_go_list = lapply(wgcna_go_files, function(f) {
	load(f)
	o = as.data.frame(go_modules)
	o$ONTOLOGY= as.character(o$ONTOLOGY)
	o$Cluster = gsub("_", "_ME", as.character(o$Cluster))
	return(o)
})

## merge
wgcna_go = do.call("rbind", wgcna_go_list)
wgcna_go$Cluster = gsub("BasoAmyg", "BLA", wgcna_go$Cluster)
wgcna_go$Cluster = gsub("MedialAmyg", "MeA", wgcna_go$Cluster)
wgcna_go$Cluster = gsub("DLPFC", "dlPFC", wgcna_go$Cluster)
wgcna_go$Cluster = gsub("Amyg", "Amygdala", wgcna_go$Cluster)

### match back
wgcna_go_bp = wgcna_go[wgcna_go$ONTOLOGY == "BP",]
wgcna_go_bp_sort = wgcna_go_bp[order(wgcna_go_bp$pvalue),]
wgcna_go_bp_top = wgcna_go_bp_sort[!duplicated(wgcna_go_bp_sort$Cluster),]
wgcna_go_bp_top = wgcna_go_bp_top[match(region_stats$Module_Name, wgcna_go_bp_top$Cluster),]
colnames(wgcna_go_bp_top) = paste0("GOBP_",colnames(wgcna_go_bp_top))

wgcna_go_mf = wgcna_go[wgcna_go$ONTOLOGY == "MF",]
wgcna_go_mf_sort = wgcna_go_mf[order(wgcna_go_mf$pvalue),]
wgcna_go_mf_top = wgcna_go_mf_sort[!duplicated(wgcna_go_mf_sort$Cluster),]
wgcna_go_mf_top = wgcna_go_mf_top[match(region_stats$Module_Name, wgcna_go_mf_top$Cluster),]
colnames(wgcna_go_mf_top) = paste0("GOMF_",colnames(wgcna_go_mf_top))

## concat gene symbols 

### join everything
region_stats = cbind(region_stats, wgcna_go_bp_top[,c(3,4,7)], wgcna_go_mf_top[,c(3,4,7)])
write.csv(region_stats, "tables/suppTable_WGCNA_enrichments.csv", row.names=FALSE)

########################
# numbers for paper ####
########################

## filters
region_stats[p.adjust(region_stats$PTSD_Pval, "fdr") < 0.05 & region_stats$PTSD_OR < 1,]
region_stats[p.adjust(region_stats$PTSD_Pval, "fdr") < 0.05 & region_stats$PTSD_OR > 1 &
	p.adjust(region_stats$PTSD_p, "fdr") < 0.05,]
region_stats[p.adjust(region_stats$MDD_Pval, "fdr") < 0.05 & region_stats$MDD_OR < 1,]
region_stats[p.adjust(region_stats$MDD_Pval, "fdr") < 0.05 & region_stats$MDD_OR > 1,]

table(p.adjust(region_stats$PTSD_Pval, "fdr") < 0.05 & region_stats$PTSD_OR > 1, 
	p.adjust(region_stats$MDD_Pval, "fdr") < 0.05 & region_stats$MDD_OR > 1)

## subset for main table
either_sig = region_stats[(p.adjust(region_stats$PTSD_Pval, "fdr") < 0.05 & region_stats$PTSD_OR > 1) | 
	(p.adjust(region_stats$MDD_Pval, "fdr") < 0.05 & region_stats$MDD_OR > 1) ,]

either_sig_pmat = either_sig[,c("Module_Name", "numGenes", 
	"PTSD_Pval", "MDD_Pval",
	"PTSD_p", "MDD_p","GOBP_Description","GOBP_pvalue" )]
write.csv(either_sig_pmat, "tables/mainTable_WGCNA_enrichment_anySig.csv",
	row.names=FALSE)

