#######
library(jaffelab)
library(limma)
library(SummarizedExperiment)
library(RColorBrewer)
library(GGally)
library(EnhancedVolcano)
dir.create("figures/figure4/")

## helper function
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
							"MeA", "dACC", "dlPFC"))				
	if(output == "wide") dd = dcast(dd, Dataset ~ Comparison)
	return(dd)
}

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

# #and counts
load("../../../count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")

## and get "cleaned" expression
load("../../Data/rdas/General/ModelMatrices/Main/AllRegions/PTSD_qsvs_Regionintxn.Rdata")
geneExprs = log2(recount::getRPKM(rse_gene, "Length") + 1)
geneExprs = geneExprs[rownames(deg_stats),]
cleanExprs = cleaningY(geneExprs, modQsva[,c(1:6, 21:26, 7:20,27:46)], P=12)

rse_gene$Group = factor(rse_gene$Group, levels = c("Control", "PTSD", "MDD"))
rse_gene$Region[rse_gene$Region == "BasoAmyg"] = "BLA"
rse_gene$Region[rse_gene$Region == "MedialAmyg"] = "MeA"
rse_gene$Region[rse_gene$Region == "DLPFC"] = "dlPFC"
rse_gene$Region = factor(rse_gene$Region, levels = c("dACC", "dlPFC", "BLA","MeA"))
subregion_pal = brewer.pal(6, "Paired")
names(subregion_pal) = levels(rse_gene$Region )

g = factor(paste0(rse_gene$Region, ":", rse_gene$Group),
	c("dACC:Control", "dACC:PTSD","dACC:MDD", 
	"dlPFC:Control", "dlPFC:PTSD", "dlPFC:MDD", 
	"BLA:Control", "BLA:PTSD", "BLA:MDD",
	"MeA:Control","MeA:PTSD", "MeA:MDD"))
	
########################
## figure 4 - mdd effects 

## get some summaries
numSig = vec_to_tab(colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.05))
colnames(numSig)[3] = "Count_FDR05"
numSig$Count_FDR10 = colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.1)
numSig[numSig$Comparison == "MDD",]

region_pal = c("darkblue", "darkgreen")
names(region_pal) = c("cortex","amygdala")

#####################
## MDD Volcanos #####
#####################

### FANCY volcano for cortex #######

pdf("figures/figure4/supp_volcano_cortex_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$Cortex_adjPVal_MDD < 0.1],
	x = 'Cortex_logFC_MDD',y = 'Cortex_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$Cortex_PValue_MDD[
		deg_stats$Cortex_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=c('darkgrey', 'darkgrey', 'darkgrey', region_pal[["cortex"]]),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

### FANCY volcano for amyg
pdf("figures/figure4/supp_volcano_amygdala_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$Amygdala_adjPVal_MDD < 0.1],
	x = 'Amygdala_logFC_MDD',y = 'Amygdala_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$Amygdala_adjPVal_MDD[
		deg_stats$Amygdala_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=c('darkgrey', 'darkgrey', 'darkgrey', region_pal[["amygdala"]]),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()


#############################
## fancy by subregion #######
pdf("figures/figure4/supp_volcano_dlpfc_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$dlPFC_adjPVal_MDD < 0.1],
	x = 'dlPFC_logFC_MDD',y = 'dlPFC_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$dlPFC_PValue_MDD[
		deg_stats$dlPFC_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["dlPFC"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure4/supp_volcano_dacc_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$dACC_adjPVal_MDD < 0.1],
	x = 'dACC_logFC_MDD',y = 'dACC_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$dACC_PValue_MDD[
		deg_stats$dACC_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["dACC"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure4/supp_volcano_mea_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$MeA_adjPVal_MDD < 0.1],
	x = 'MeA_logFC_MDD',y = 'MeA_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$MeA_PValue_MDD[
		deg_stats$MeA_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["MeA"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure4/supp_volcano_bla_mdd_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$BLA_adjPVal_MDD < 0.1],
	x = 'BLA_logFC_MDD',y = 'BLA_PValue_MDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$BLA_PValue_MDD[
		deg_stats$BLA_adjPVal_MDD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["BLA"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

##############
## scatters supp
######

## t
mdd_ts = deg_stats[,grep("_t_MDD$", colnames(deg_stats))]
mdd_ts = mdd_ts[,c(1,6,7,2,4,5,3)]
colnames(mdd_ts) = ss(colnames(mdd_ts), "_")

cc_mdd = cor(mdd_ts)
round(cc_mdd,3)

pdf("figures/figure4/supp_ggpairs_tstats.pdf",h=6,w=6)
ggpairs(mdd_ts,
	lower = list(continuous = wrap("points", size=0.2)),
	upper = list(continuous = wrap('cor', size = 5, col = "black"))) + 
	theme_bw() + theme(axis.text = element_text(size = 8))
dev.off()

## log fc
mdd_lfcs = deg_stats[,grep("_logFC_MDD$", colnames(deg_stats))]
mdd_lfcs = mdd_lfcs[,c(1,6,7,2,4,5,3)]
colnames(mdd_lfcs) = ss(colnames(mdd_lfcs), "_")

cc_lfcs = cor(mdd_lfcs)
round(cc_lfcs,3)

pdf("figures/figure4/supp_ggpairs_lfcs.pdf",h=6,w=6)
ggpairs(mdd_lfcs,
	lower = list(continuous = wrap("points", size=0.2)),
	upper = list(continuous = wrap('cor', size = 5, col = "black"))) + 
	theme_bw() + theme(axis.text = element_text(size = 8))
dev.off()

#######################
#### ptsd vs mdd ######
#######################

#### correlations
ptsd_ts = deg_stats[,grep("_t_PTSD$", colnames(deg_stats))]
ptsd_ts = ptsd_ts[,c(1,6,7,2,4,5,3)]
colnames(ptsd_ts) = ss(colnames(ptsd_ts), "_")

mdd_ptsd_corr= cor(mdd_ts, ptsd_ts)
cors = signif(diag(mdd_ptsd_corr),3)
names(cors) = colnames(ptsd_ts)
cors
range(cors)

## show cortex for main panel, better labels
cols = rep("lightgrey", nrow(deg_stats))
cols[deg_stats$Cortex_PValue_MDD < 0.005 & deg_stats$Cortex_PValue_PTSD > 0.005] = "blue"
cols[deg_stats$Cortex_PValue_MDD > 0.005 & deg_stats$Cortex_PValue_PTSD < 0.005] = "red"
cols[deg_stats$Cortex_PValue_MDD < 0.005 & deg_stats$Cortex_PValue_PTSD < 0.005] = "purple"

pdf("figures/figure4/A_cortex_scatter_mdd_vs_ptsd.pdf")
par(mar=c(5,6,2,2),cex.axis=2.3,cex.lab=2.3)
plot(deg_stats$Cortex_t_PTSD, deg_stats$Cortex_t_MDD,
	pch = 21, bg=cols,
	xlab = "PTSD Effect", ylab = "MDD Effect",
	xlim = c(-5.2,5.2), ylim = c(-5.2,5.2))
text(x=0,y=4.75,"Cortex",cex=3)
abline(h=c(-2.82, 2.82), v=c(-2.82, 2.82), lty=2,col="blue")
abline(0,1, lty=2,col="black")
dev.off()

## scatterplots for supp
pdf("figures/figure4/scatterplots_mdd_vs_ptsd.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
for(i in 1:ncol(mdd_ts)) {
	plot(ptsd_ts[,i], mdd_ts[,i], pch = 21, bg="grey",
		xlab = "PTSD vs Control (T-stat)", 
		ylab = "MDD vs Control (T-stat)",
		xlim = c(-6,6), ylim = c(-6,6))
	legend("topleft", colnames(ptsd_ts)[i],cex=3,bty="n")
	abline(h=c(-2.8, 2.8), v=c(-2.8, 2.8), lty=2,col="blue")
	abline(0,1, lty=2,col="blue")
}
dev.off()


##### overlap among significant genes
geneSets_PTSD = with(deg_stats, 
	list(Cortex_either = Cortex_PValue_PTSD < 0.005,
		dlPFC_either = dlPFC_PValue_PTSD < 0.005,
		dACC_either = dACC_PValue_PTSD < 0.005,
		Amygdala_either = Amygdala_PValue_PTSD < 0.005,
		BLA_either = BLA_PValue_PTSD < 0.005,
		MeA_either = MeA_PValue_PTSD < 0.005))
geneSets_MDD = with(deg_stats, 
	list(Cortex_either = Cortex_PValue_MDD < 0.005,
		dlPFC_either = dlPFC_PValue_MDD < 0.005,
		dACC_either = dACC_PValue_MDD < 0.005,
		Amygdala_either = Amygdala_PValue_MDD < 0.005,
		BLA_either = BLA_PValue_MDD < 0.005,
		MeA_either = MeA_PValue_MDD < 0.005))
	
mdd_ptsd_overlap_tables = mapply(function(mdd, ptsd) {
	tt = table(mdd, ptsd, dnn = c("MDD", "PTSD"))
}, geneSets_MDD, geneSets_PTSD,SIMPLIFY=FALSE)

## ptsd dir
sapply(mdd_ptsd_overlap_tables, function(x) {
	prop.table(x, 2)[2,2]
})
## mdd dir
sapply(mdd_ptsd_overlap_tables, function(x) {
	prop.table(x, 1)[2,2]
})
fisher_list = lapply(mdd_ptsd_overlap_tables, fisher.test)
sapply(fisher_list, "[[", "p.value")

#################
## MDD vs PTSD ##
#################
numSig$P005= colSums(deg_stats[,grep("PValue", colnames(deg_stats))] < 0.005)

numSig[numSig$Comparison == "PTSDvsMDD",]

deg_stats[deg_stats$Cortex_adjPVal_PTSDvsMDD < 0.1,]
deg_stats[deg_stats$MeA_adjPVal_PTSDvsMDD < 0.1,]
deg_stats[deg_stats$dlPFC_adjPVal_PTSDvsMDD < 0.1,]

### FANCY volcano for cortex
pdf("figures/figure4/B_volcano_cortex_PTSDvsMDD_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding") &
		deg_stats$Cortex_PValue_PTSDvsMDD < 0.005],
	x = 'Cortex_logFC_PTSDvsMDD',y = 'Cortex_PValue_PTSDvsMDD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = 0.005,
	xlim = c(-1.15,1.15),ylim = c(0,5.25),
	FCcutoff = 0, title = "", subtitle = "",
	col=c('darkgrey', 'darkgrey', 'darkgrey', region_pal[["cortex"]]),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

################
## PTSDonly ####
################

cols = rep("lightgrey", nrow(deg_stats))
cols[deg_stats$Cortex_PValue_onlyPTSD < 0.005 & deg_stats$Cortex_PValue_PTSD > 0.005] = "blue"
cols[deg_stats$Cortex_PValue_onlyPTSD > 0.005 & deg_stats$Cortex_PValue_PTSD < 0.005] = "red"
cols[deg_stats$Cortex_PValue_onlyPTSD < 0.005 & deg_stats$Cortex_PValue_PTSD < 0.005] = "purple"


pdf("figures/figure4/D_cortex_scatter_ptsd_vs_only.pdf")
par(mar=c(5,6,2,2),cex.axis=2.3,cex.lab=2.3)
plot(deg_stats$Cortex_t_PTSD, deg_stats$Cortex_t_onlyPTSD,
	pch = 21, bg=cols,
	xlab = "PTSD-Control Effect", ylab = "PTSD-Specific Effect",
	xlim = c(-5.2,5.2), ylim = c(-5.2,5.2))
text(x=0,y=4.75,"Cortex",cex=3)
abline(h=c(-2.82, 2.82), v=c(-2.82, 2.82), lty=2,col=region_pal[["cortex"]])
abline(0,1, lty=2,col="black")
dev.off()
