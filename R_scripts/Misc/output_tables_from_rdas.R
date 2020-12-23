##
library(GenomicRanges)

############
## GENE ####
############

deg_files = c(Cortex = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",
	Amygdala = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",
	Joint = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda",
	BLA= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda",
	MeA= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda",
	dACC= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda",
	DLPFC= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")

deg_names = sapply(deg_files, load, verbose=TRUE, envir = .GlobalEnv)
deg_list = mapply(function(f, n) {
	load(f)
	return(get(n))
}, f = deg_files, n = deg_names)

sapply(deg_list, dim)

## get the stats
deg_stat_list = lapply(deg_list, function(x) {
	ind = sort(c(grep("Val", colnames(x)),grep("^t", colnames(x)),
				grep("^log", colnames(x))))

	x[,ind]
})	
## merge
deg_stats = do.call("cbind", deg_stat_list)
## drop regoin associations
deg_stats = deg_stats[,!grepl("_Region", colnames(deg_stats))]

## rename
colnames(deg_stats) = gsub("P.Value", "PValue", colnames(deg_stats) )
colnames(deg_stats) = gsub("adj.P.Val", "adjPVal", colnames(deg_stats) )
colnames(deg_stats) = gsub("\\.", "_", colnames(deg_stats) )

## add annotation
map = deg_list[[1]][,c(165, 163, 164, 161,166)]

deg_stats = cbind(map, deg_stats)
write.csv(deg_stats, file = gzfile("../csvs/merged_deg_stats_allComparisons.csv.gz"))

############
## EXON ####
############

dee_files = c(Cortex = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda",
	Amygdala = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda",
	Joint = "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/AllRegions/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda",
	BLA= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/BasoAmyg/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda",
	MeA= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/MedialAmyg/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda",
	dACC= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/dACC/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda",
	DLPFC= "/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/DEG/DLPFC/exonStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda")

dee_names = sapply(dee_files, load, verbose=TRUE, envir = .GlobalEnv)
dee_list = mapply(function(f, n) {
	load(f)
	return(get(n))
}, f = dee_files, n = dee_names)

sapply(dee_list, dim)

## get the stats
dee_stat_list = lapply(dee_list, function(x) {
	ind = sort(c(grep("Val", colnames(x)),grep("^t", colnames(x)),
				grep("^log", colnames(x))))

	x[,ind]
})	
## merge
dee_stats = do.call("cbind", dee_stat_list)
## drop regoin associations
dee_stats = dee_stats[,!grepl("_Region", colnames(dee_stats))]

## rename
colnames(dee_stats) = gsub("P.Value", "PValue", colnames(dee_stats) )
colnames(dee_stats) = gsub("adj.P.Val", "adjPVal", colnames(dee_stats) )
colnames(dee_stats) = gsub("\\.", "_", colnames(dee_stats) )

## add annotation
map = dee_list[[1]][,c(212:225)]

dee_stats = cbind(map, dee_stats)
save(dee_stats, file = "../../Misc_large_files/merged_dee_stats_allComparisons.Rdata")