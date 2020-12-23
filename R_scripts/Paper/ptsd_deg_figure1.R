#######
library(jaffelab)
library(limma)
library(SummarizedExperiment)
library(RColorBrewer)
library(GGally)
library(EnhancedVolcano)

## create folder
dir.create("figures")
dir.create("tables")

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
	
## samples
table(rse_gene$Group[!duplicated(rse_gene$BrNum)])
# Control     MDD    PTSD
    # 109     109     107

########################
## figure 1 - ptsd effects 

## get some summaries
numSig = vec_to_tab(colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.05))
colnames(numSig)[3] = "Count_FDR05"
numSig$Count_FDR10 = colSums(deg_stats[,grep("adjP", colnames(deg_stats))] < 0.1)
numSig[numSig$Comparison == "PTSD",]

table(deg_stats$dlPFC_adjPVal_PTSD < 0.05 | deg_stats$dACC_adjPVal_PTSD < 0.05 )
table(deg_stats$dlPFC_adjPVal_PTSD < 0.1 | deg_stats$dACC_adjPVal_PTSD < 0.1 )
deg_stats$Symbol[which(deg_stats$dlPFC_adjPVal_PTSD < 0.1 & deg_stats$dACC_adjPVal_PTSD < 0.1 )]
table(deg_stats$BLA_adjPVal_PTSD < 0.1 | deg_stats$MeA_adjPVal_PTSD < 0.1 )
deg_stats$Symbol[which(deg_stats$BLA_adjPVal_PTSD < 0.1 & deg_stats$MeA_adjPVal_PTSD < 0.1 )]

## interaction effects
deg_stats_cortex = deg_stats[deg_stats$Cortex_adjPVal_PTSD < 0.1,
	c(1:6, grep("Cortex",colnames(deg_stats)) )]
min(deg_stats_cortex$Cortex_PValue_rintxnPTSD)
cat(deg_stats_cortex$Symbol[deg_stats_cortex$Cortex_PValue_rintxnPTSD < 0.01],sep=", ")
deg_stats_amygdala = deg_stats[deg_stats$Amygdala_adjPVal_PTSD < 0.1,
	c(1:6, grep("Amygdala",colnames(deg_stats)) )]
min(deg_stats_amygdala$Amygdala_PValue_rintxnPTSD)
cat(deg_stats_amygdala$Symbol[deg_stats_cortex$Amygdala_PValue_rintxnPTSD < 0.01],sep=", ")

region_pal = c("darkblue", "darkgreen")
names(region_pal) = c("cortex","amygdala")

dir.create("figures/figure1/")

### FANCY?
pdf("figures/figure1/A_volcano_cortex_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$Cortex_adjPVal_PTSD < 0.05],
	x = 'Cortex_logFC_PTSD',y = 'Cortex_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$Cortex_PValue_PTSD[
		deg_stats$Cortex_adjPVal_PTSD < 0.05]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=c('darkgrey', 'darkgrey', 'darkgrey', region_pal[["cortex"]]),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure1/B_volcano_amygdala_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$Amygdala_adjPVal_PTSD < 0.05],
	x = 'Amygdala_logFC_PTSD',y = 'Amygdala_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = 5e-6,
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=c('darkgrey', 'darkgrey', 'darkgrey', region_pal[["amygdala"]]),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

#############################
## fancy by subregion #######
pdf("figures/figure1/supp_volcano_dlpfc_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$dlPFC_adjPVal_PTSD < 0.1],
	x = 'dlPFC_logFC_PTSD',y = 'dlPFC_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$dlPFC_PValue_PTSD[
		deg_stats$dlPFC_adjPVal_PTSD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["dlPFC"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure1/supp_volcano_dacc_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$dACC_adjPVal_PTSD < 0.1],
	x = 'dACC_logFC_PTSD',y = 'dACC_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$dACC_PValue_PTSD[
		deg_stats$dACC_adjPVal_PTSD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["dACC"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure1/supp_volcano_mea_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$MeA_adjPVal_PTSD < 0.1],
	x = 'MeA_logFC_PTSD',y = 'MeA_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$MeA_PValue_PTSD[
		deg_stats$MeA_adjPVal_PTSD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["MeA"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure1/supp_volcano_bla_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$BLA_adjPVal_PTSD < 0.1],
	x = 'BLA_logFC_PTSD',y = 'BLA_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$BLA_PValue_PTSD[
		deg_stats$BLA_adjPVal_PTSD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "", subtitle = "",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', subregion_pal["BLA"])),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()

pdf("figures/figure1/supp_volcano_joint_fancy.pdf", useDingbats= FALSE)
EnhancedVolcano(deg_stats,
	lab = deg_stats$Symbol,
	selectLab = deg_stats$Symbol[
		deg_stats$gene_type %in% c("protein_coding","lincRNA") &
		deg_stats$Joint_adjPVal_PTSD < 0.1],
	x = 'Joint_logFC_PTSD',y = 'Joint_PValue_PTSD',
	drawConnectors = FALSE,
	legendPosition = 'none',  labSize = 4.0,
	pCutoff = max(deg_stats$Joint_PValue_PTSD[
		deg_stats$Joint_adjPVal_PTSD < 0.1]),
	xlim = c(-1.15,1.15),ylim = c(0,8),
	FCcutoff = 0, title = "Joint Data", subtitle = "PTSD vs Control Effects",
	col=as.character(c('darkgrey', 'darkgrey', 'darkgrey', "black")),
	caption = paste0('Total = ', nrow(deg_stats), ' genes'))
dev.off()


###########################
## single gene examples ##
##########################

pdf("figures/figure1/C_gene_boxplots_cleaned.pdf",w=6,h=5)
gene_to_plot = which(deg_stats$Amygdala_adjPVal_PTSD < 0.1 | 
	deg_stats$Cortex_adjPVal_PTSD < 0.1)
palette(subregion_pal)
par(mar=c(7,6,2,2),cex.axis=1.8,cex.lab=1.8,cex.main= 2)
for( i in gene_to_plot) {
	boxplot(cleanExprs[i,] ~ g, outline=FALSE, xlab="",las=3,
		ylim = quantile(cleanExprs[i,],c(0,0.995)),
		names= ss(levels(g), ":", 2),
		main = deg_stats$Symbol[i],ylab = "Adjusted Exprs (log2)")
	points(cleanExprs[i,] ~ jitter(as.numeric(g),amount=0.15),
		pch = 20+as.numeric(rse_gene$Group),
		bg = rse_gene$Region)
	text(x = c(2,5,8,11), y=min(cleanExprs[i,]), 
		levels(rse_gene$Region),cex=2)
	abline(v=c(3.5,6.5,9.5), lty=1.5)
}
dev.off()
	
pdf("figures/figure1/C_gene_boxplots_obs.pdf",w=6,h=5)
gene_to_plot = which(deg_stats$Amygdala_adjPVal_PTSD < 0.1 | 
	deg_stats$Cortex_adjPVal_PTSD < 0.1)
palette(subregion_pal)
par(mar=c(7,6,2,2),cex.axis=1.8,cex.lab=1.8,cex.main= 2)
for( i in gene_to_plot) {
	boxplot(geneExprs[i,] ~ g, outline=FALSE, xlab="",las=3,
		ylim = quantile(geneExprs[i,],c(0,0.995)),
		names= ss(levels(g), ":", 2),
		main = deg_stats$Symbol[i],ylab = "log2[RPKM+1]")
	points(geneExprs[i,] ~ jitter(as.numeric(g),amount=0.15),
		pch = 20+as.numeric(rse_gene$Group),
		bg = rse_gene$Region)
	text(x = c(2,5,8,11), y=min(geneExprs[i,]), 
		levels(rse_gene$Region),cex=1.5)
	abline(v=c(3.5,6.5,9.5), lty=2)
}
dev.off()

##############
## scatters supp
######

## t
ptsd_ts = deg_stats[,grep("_t_PTSD$", colnames(deg_stats))]
ptsd_ts = ptsd_ts[,c(1,6,7,2,4,5,3)]
colnames(ptsd_ts) = ss(colnames(ptsd_ts), "_")

cc_ptsd = cor(ptsd_ts)
round(cc_ptsd,3)

pdf("figures/figure1/supp_ggpairs_tstats.pdf",h=6,w=6)
ggpairs(ptsd_ts,
	lower = list(continuous = wrap("points", size=0.2)),
	upper = list(continuous = wrap('cor', size = 5, col = "black"))) + 
	theme_bw() + theme(axis.text = element_text(size = 8))
dev.off()

## log fc
ptsd_lfcs = deg_stats[,grep("_logFC_PTSD$", colnames(deg_stats))]
ptsd_lfcs = ptsd_lfcs[,c(1,6,7,2,4,5,3)]
colnames(ptsd_lfcs) = ss(colnames(ptsd_lfcs), "_")

cc_lfcs = cor(ptsd_lfcs)
round(cc_lfcs,3)

pdf("figures/figure1/supp_ggpairs_lfcs.pdf",h=6,w=6)
ggpairs(ptsd_lfcs,
	lower = list(continuous = wrap("points", size=0.2)),
	upper = list(continuous = wrap('cor', size = 5, col = "black"))) + 
	theme_bw() + theme(axis.text = element_text(size = 8))
dev.off()

#######################
## sensitivity stuff ##
#######################

sex_effect_csvs = list(dACC_F = "../../csvs/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_females_threeGroup.csv.gz",
						dACC_M = "../../csvs/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_males_threeGroup.csv.gz",
						dlPFC_F = "../../csvs/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_females_threeGroup.csv.gz",
						dlPFC_M = "../../csvs/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_males_threeGroup.csv.gz",
						BLA_F = "../../csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_females_threeGroup.csv.gz",
						BLA_M = "../../csvs/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_males_threeGroup.csv.gz",
						MeA_F = "../../csvs/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_females_threeGroup.csv.gz",
						MeA_M = "../../csvs/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_males_threeGroup.csv.gz")
sex_effects_list = lapply(sex_effect_csvs, read.csv,as.is=TRUE,row.names=1)						
ptsd_sex_effect_list = lapply(sex_effects_list, function(x) {
	x[rownames(deg_stats),grep("_PTSD$", colnames(x))][,c(11,13,14,15)]
})
	
ptsd_sex_ts = sapply(ptsd_sex_effect_list, "[[", "t_PTSD")
ptsd_sex_pvals = sapply(ptsd_sex_effect_list, "[[", "P.Value_PTSD")
ptsd_sex_fdrs = sapply(ptsd_sex_effect_list, "[[", "adj.P.Val_PTSD")
rownames(ptsd_sex_ts) = rownames(ptsd_sex_pvals) = rownames(ptsd_sex_fdrs) = rownames(deg_stats)

## check all
cor(ptsd_ts, ptsd_sex_ts)
signif(cor(ptsd_ts[,names(subregion_pal)], ptsd_sex_ts),2)
signif(cor(ptsd_ts[rownames(deg_stats_cortex),names(subregion_pal)], 
	ptsd_sex_ts[rownames(deg_stats_cortex),]),2)

## cortex sex effects
pdf("figures/figure1/supp_sex_ptsd_diffs.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(subregion_pal)
plot(ptsd_sex_ts[rownames(deg_stats_cortex),c("dACC_F", "dACC_M")],
	pch = 21, bg = 1, cex= -log10(deg_stats_cortex$Cortex_PValue_PTSD)/2,
	xlab = "Female PTSD Effect", ylab = "Male PTSD Effect",
	ylim = c(-5,5),xlim = c(-5,5))
text(x=-3,y=3, "dACC",cex=3)
abline(0,1,lty=1,lwd=2)
abline(h=0,v=0,lty=2)

plot(ptsd_sex_ts[rownames(deg_stats_cortex),c("dlPFC_F", "dlPFC_M")],
	pch = 21, bg = 2, cex= -log10(deg_stats_cortex$Cortex_PValue_PTSD)/2,
	xlab = "Female PTSD Effect", ylab = "Male PTSD Effect",
	ylim = c(-5,5),xlim = c(-5,5))
text(x=-3,y=3, "dlPFC",cex=3)
abline(0,1,lty=1,lwd=2)
abline(h=0,v=0,lty=2)
dev.off()

fvmt = sapply(split0(ss(colnames(ptsd_sex_ts),"_")), function(x) {
	ptsd_sex_ts[,x[1]] - ptsd_sex_ts[,x[2]]
})
colnames(fvmt) = paste0(colnames(fvmt), "_FvM") 
## for fancy plots
sex_checks = as.data.frame(cbind(ptsd_ts, fvmt))
sex_checks$Symbol = deg_stats$Symbol
sex_checks$gene_type = deg_stats$gene_type
sex_checks_cortex = sex_checks[rownames(deg_stats_cortex),]


## cortex sex effects
pdf("figures/figure1/supp_sex_ptsd_mastyle.pdf")
set.seed(3)
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(subregion_pal)
plot(sex_checks_cortex$dACC, sex_checks_cortex$dACC_FvM,
	pch = 21, bg = 1, 
	xlab = "Overall PTSD Effect", ylab = "Male    Sex Bias     Female",
	ylim = c(-3,3),xlim = c(-6,6))
ind = which(abs( sex_checks_cortex$dACC_FvM) > 2 & 
	sex_checks_cortex$gene_type %in% c("protein_coding","lincRNA"))
text(sex_checks_cortex$dACC[ind], sex_checks_cortex$dACC_FvM[ind],
	sex_checks_cortex$Symbol[ind],cex=0.6,
	pos=sample(1:4, length(ind),replace=TRUE))
	
plot(sex_checks_cortex$dlPFC, sex_checks_cortex$dlPFC_FvM,
	pch = 21, bg = 2, 
	xlab = "Overall PTSD Effect", ylab = "Male    Sex Bias    Female",
	ylim = c(-3,3),xlim = c(-6,6))
ind = which(abs( sex_checks_cortex$dlPFC_FvM) > 2 & 
	sex_checks_cortex$gene_type %in% c("protein_coding","lincRNA"))
text(sex_checks_cortex$dlPFC[ind], sex_checks_cortex$dlPFC_FvM[ind],
	sex_checks_cortex$Symbol[ind],cex=0.6,
	pos=sample(1:4, length(ind),replace=TRUE))

dev.off()


#####################
##  check exons #####
#####################

load("../../Misc_large_files/merged_dee_stats_allComparisons.Rdata")
colnames(dee_stats) = gsub("DLPFC", "dlPFC", colnames(dee_stats))
dee_stats = as.data.frame(dee_stats)

## get some summaries
numSig_exons = vec_to_tab(colSums(dee_stats[,grep("adjP", colnames(dee_stats))] < 0.05))
colnames(numSig_exons)[3] = "Count_FDR05"
numSig_exons$Count_FDR10 = colSums(dee_stats[,grep("adjP", colnames(dee_stats))] < 0.1)
numSig_exons[numSig_exons$Comparison == "PTSD",]

numSig$Count_FDR05_Exon = numSig_exons$Count_FDR05
numSig$Count_FDR10_Exon = numSig_exons$Count_FDR10
numSig[numSig$Comparison == "PTSD",]

table(dee_stats$Cortex_adjPVal_PTSD < 0.05 | dee_stats$Amygdala_adjPVal_PTSD < 0.05 )
unique(dee_stats$Symbol[dee_stats$Cortex_adjPVal_PTSD < 0.05 | dee_stats$Amygdala_adjPVal_PTSD < 0.05 ])
xx = dee_stats[dee_stats$Cortex_adjPVal_PTSD < 0.05 | dee_stats$Amygdala_adjPVal_PTSD < 0.05 ,]
xx[order(xx$Cortex_adjPVal_PTSD),1:18 ]

table(dee_stats$dlPFC_adjPVal_PTSD < 0.1 | dee_stats$dACC_adjPVal_PTSD < 0.1 )
unique(dee_stats$Symbol[which(dee_stats$dlPFC_adjPVal_PTSD < 0.1 | dee_stats$dACC_adjPVal_PTSD < 0.1 )])
unique(dee_stats$Symbol[which(dee_stats$dlPFC_adjPVal_PTSD < 0.05 | dee_stats$dACC_adjPVal_PTSD < 0.05 )])
table(dee_stats$BLA_adjPVal_PTSD < 0.1 | dee_stats$MeA_adjPVal_PTSD < 0.1 )
dee_stats$Symbol[which(dee_stats$BLA_adjPVal_PTSD < 0.1 & dee_stats$MeA_adjPVal_PTSD < 0.1 )]
