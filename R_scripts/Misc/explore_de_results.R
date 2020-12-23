###
library(reshape2)
library(jaffelab)
library(limma)
library(qvalue)

vec_to_tab = function(vec, output = "long") {
	dd = data.frame(Count = vec)
	dd$Dataset = ss(rownames(dd), "_", 1)
	dd$Comparison = ss(rownames(dd), "_", 3)
	rownames(dd) = NULL
	dd = dd[,c("Dataset", "Comparison", "Count")]
	dd$Comparison = factor(dd$Comparison,	
			levels = c("PTSD", "MDD", "ANOVA", "onlyPTSD", 
				"PTSDvsMDD", "rintxnMDD", "rintxnPTSD", "rintxnonlyPTSD"))
	dd$Dataset = factor(dd$Dataset,	
			levels = c("Cortex", "Amygdala", "Joint", "BLA",
							"MeA", "dACC", "DLPFC"))				
	if(output == "wide") dd = dcast(dd, Dataset ~ Comparison)
	return(dd)
}

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)

## get some summaries
numSig = vec_to_tab(colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.05))
colnames(numSig)[3] = "Count_FDR05"
numSig$Count_FDR10 = colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.1)

vec_to_tab(colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.05), "wide" )
vec_to_tab(colSums(deg_stats[,grep("PValue", colnames(deg_stats))] < 0.01), "wide" )

## genes in each specific analysis
ptsd_fdrs = deg_stats[,grep("adjPVal_PTSD$", colnames(deg_stats))]
sigIndex_ptsd = apply(ptsd_fdrs, 2, function(x) rownames(ptsd_fdrs)[x < 0.1])
## cross model checks

## with q-value?
qval = apply(deg_stats[,grep("PValue", colnames(deg_stats))], 2, 
	function(x) qvalue(x)$qvalues)
numSig$Count_StoreyFDR05 = colSums(qval < 0.05)
numSig$Count_StoreyFDR10 = colSums(qval < 0.1)
#######################################
### check replication with PTSD #######
yale = read.csv("../../../Yale_SuppTables_cleaned/S1-SupplementaryTable-PTSD-DEGs-fixedgenenames.csv",
	row.names = 1, as.is=TRUE)

## overall
mm_yale = match(deg_stats$ensemblID, yale$Geneid)
cc_yale = cor(deg_stats[!is.na(mm_yale),grep("logFC_PTSD$", colnames(deg_stats))],
	yale[mm_yale[!is.na(mm_yale)], grep("log2", colnames(yale))],
		use = "pair")
colnames(cc_yale) = ss(colnames(cc_yale), "\\.", 2)
rownames(cc_yale) = ss(rownames(cc_yale), "_")
signif(cc_yale,3)

###########################
## libd significant gene ##
g_union_libd = deg_stats$ensemblID[rowSums(deg_stats[,
	grep("adjPVal_PTSD$", colnames(deg_stats))] < 0.1) > 0]
mm_yale_sig = match(g_union_libd, yale$Geneid)
cc_yale_sig = cor(deg_stats[match(g_union_libd,deg_stats$ensemblID),
			grep("logFC_PTSD$", colnames(deg_stats))],
	yale[mm_yale_sig, grep("log2", colnames(yale))],
		use = "pair")
colnames(cc_yale_sig) = ss(colnames(cc_yale_sig), "\\.", 2)
rownames(cc_yale_sig) = ss(rownames(cc_yale_sig), "_")
signif(cc_yale_sig,3)

############################
## yale significant genes ##
g_union_yale = yale$Geneid[rowSums(yale[,
	grep("padj", colnames(yale))] < 0.1, na.rm=TRUE) > 0]
mm_libd_sig = match(g_union_yale, deg_stats$ensemblID)
cc_libd_sig = cor(deg_stats[mm_libd_sig,
			grep("logFC_PTSD$", colnames(deg_stats))],
	yale[match(g_union_yale,yale$Geneid), grep("log2", colnames(yale))],
		use = "pair")
colnames(cc_libd_sig) = ss(colnames(cc_libd_sig), "\\.", 2)
rownames(cc_libd_sig) = ss(rownames(cc_libd_sig), "_")
signif(cc_libd_sig,3)

## which approach had higher replication?
cc_yale_sig > cc_libd_sig
mean(cc_yale_sig > cc_libd_sig) # 0.85

#####################################
### read in mdd replication stats ###

load("/dcl02/lieber/ajaffe/PublicData/MDD_Nestler_New/mdd_effects_3models_topTables.Rdata")

mdd_t_nestler_adj = sapply(topTableList_other, function(x) x[rownames(deg_stats),"t"])
rownames(mdd_t_nestler_adj) = rownames(deg_stats)
mdd_t_nestler_qsva = sapply(topTableList_qSVA, function(x) x[rownames(deg_stats),"t"])
rownames(mdd_t_nestler_qsva) = rownames(deg_stats)

## overall, qSVA
signif(cor(deg_stats[,grep("t_MDD$", colnames(deg_stats))], mdd_t_nestler_qsva),3)

## significant in ours
g_union_mdd = rowSums(deg_stats[,grep("adjPVal_MDD$", colnames(deg_stats))] < 0.1) > 0
signif(cor(deg_stats[g_union_mdd,grep("t_MDD$", colnames(deg_stats))], 
	mdd_t_nestler_qsva[g_union_mdd,]),3)

## no latent variables
signif(cor(deg_stats[,grep("t_MDD$", colnames(deg_stats))], mdd_t_nestler_adj),3)
signif(cor(deg_stats[g_union_mdd,grep("t_MDD$", colnames(deg_stats))], 
	mdd_t_nestler_adj[g_union_mdd,]),3)
### no overlap either way

####################
## ptsd vs mdd #####
####################

sig_index_any = rowSums(deg_stats[,grep("adjPVal_PTSD$", colnames(deg_stats))] < 0.1) > 0
	
cor(deg_stats[sig_index_any,grep("t_onlyPTSD",colnames(deg_stats))],
	deg_stats[sig_index_any,grep("t_MDD",colnames(deg_stats))])
cor(deg_stats[sig_index_any,grep("t_onlyPTSD",colnames(deg_stats))],
	deg_stats[sig_index_any,grep("t_PTSDvsMDD",colnames(deg_stats))])
ggpairs(deg_stats[sig_index_any,grep("t_onlyPTSD",colnames(deg_stats))],
	deg_stats[sig_index_any,grep("t_PTSDvsMDD",colnames(deg_stats))])
grep("t_MDD$",colnames(deg_stats)) 
grep("t_PTSDvsMDD",colnames(deg_stats)) 
