###

#Read in our data
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('clusterProfiler')
library(sva)
library('readxl')
library(recount)
library(limma)
library(edgeR)
library(RColorBrewer)

dir.create("figures/figure5")

## read back in merged stats
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))

# #and counts
load("../../../count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")

## and get model to "clean" expression
load("../../Data/rdas/General/ModelMatrices/Main/AllRegions/PTSD_qsvs_Regionintxn.Rdata")
geneExprs = log2(recount::getRPKM(rse_gene, "Length") + 1)
geneExprs = geneExprs[rownames(deg_stats),]

## and get model 
load("../../Data/rdas/General/ModelMatrices/Main/AllRegions/PTSD_qsvs_Regionintxn.Rdata")
## no region or interaction
modQsva_noRegion = modQsva[,c(1:3, 7:20, 27:46)]

#####################
## group labels #####
rse_gene$Group = factor(rse_gene$Group, levels = c("Control", "PTSD", "MDD"))
rse_gene$Region[rse_gene$Region == "BasoAmyg"] = "BLA"
rse_gene$Region[rse_gene$Region == "MedialAmyg"] = "MeA"
rse_gene$Region[rse_gene$Region == "DLPFC"] = "dlPFC"
rse_gene$Region = factor(rse_gene$Region, levels = c("dACC", "dlPFC", "BLA","MeA"))
subregion_pal = brewer.pal(6, "Paired")[3:6]
names(subregion_pal) = levels(rse_gene$Region )

g = factor(paste0(rse_gene$Region, ":", rse_gene$Group),
	c("dACC:Control", "dACC:PTSD","dACC:MDD", 
	"dlPFC:Control", "dlPFC:PTSD", "dlPFC:MDD", 
	"BLA:Control", "BLA:PTSD", "BLA:MDD",
	"MeA:Control","MeA:PTSD", "MeA:MDD"))
	
	
##################################
## drop to common set of brains ##
##################################

brTab = table(rse_gene$BrNum, rse_gene$Region)
all_four = rowSums(brTab == 1) == 4 
br_four = names(which(all_four))

## filter
keepIndex = rse_gene$BrNum %in% br_four
rse_gene_all4 = rse_gene[,keepIndex]
geneExprs_all4 = geneExprs[,keepIndex]
modQsva_all4 = modQsva_noRegion[keepIndex,]

## put in brnum order
sample_order = order(rse_gene_all4$BrNum, rse_gene_all4$Region)
rse_gene_all4 = rse_gene_all4[,sample_order]
geneExprs_all4 = geneExprs_all4[,sample_order]
modQsva_all4 = modQsva_all4[sample_order,]

## put br
colnames(geneExprs_all4) = rse_gene_all4$BrNum

##################################
region_pairs = as.data.frame(combn(levels(rse_gene_all4$Region), 2),stringsAsFactors=FALSE)
names(region_pairs) = apply(region_pairs,2,paste,collapse="_")
rIndexes = splitit(rse_gene_all4$Region)

### observed expression
region_coherence = sapply(region_pairs, function(r) {
	r_ind1 = rIndexes[[r[1]]]
	r_ind2 = rIndexes[[r[2]]]
	g1 = geneExprs_all4[,r_ind1]
	g2 = geneExprs_all4[,r_ind2]
	message(identical(colnames(g1),colnames(g2)))
	cc = cor(g1,g2)
	diag(cc)
})

## associations
pd = colData(rse_gene_all4)
pd = pd[!duplicated(pd$BrNum),c(2,5,6,7,8,9)]
rownames(pd) = pd$BrNum
pd = pd[rownames(region_coherence),]

coherence_effects = lapply(as.data.frame(region_coherence), function(y) {
	summary(lm(y~ Group + AgeDeath + Sex + Race, data=pd))
})
ptsd_effect = t(sapply(coherence_effects, function(x) x$coef[2,]))
mdd_effect = t(sapply(coherence_effects, function(x) x$coef[3,]))

pdf("figures/figure5/coherence_boxplots.pdf",h=6,w=6)
par(mar=c(3,6,2,2),cex.axis=2,cex.lab=2.5)
for(i in 1:ncol(region_coherence)) {
	boxplot(region_coherence[,i] ~ pd$Group,
		ylim = c(0.85,1),outline=FALSE,xlab="",
		ylab= "Coherence")
	text(x = 2, y = 1, gsub("_", "-", colnames(region_coherence)[i]),cex=2.5)
	points(region_coherence[,i] ~ jitter(as.numeric(pd$Group),amount=0.15),
		pch = 21 ,bg="grey")
}
dev.off()

### cleaned expression
region_coherence_clean = sapply(region_pairs, function(r) {
	r_ind1 = rIndexes[[r[1]]]
	r_ind2 = rIndexes[[r[2]]]
	g1 = geneExprs_all4[,r_ind1]
	g2 = geneExprs_all4[,r_ind2]
	m1 = modQsva_all4[r_ind1,]
	m2 = modQsva_all4[r_ind2,]
	## dont protect dx, but dont fit dx either
	cleang1 = cleaningY(g1,m1,P=5)
	cleang2 = cleaningY(g2,m2,P=5)
	
	message(identical(colnames(cleang1),colnames(cleang2)))
	cc = cor(cleang1,cleang2)
	diag(cc)
})
coherence_effects_clean = lapply(as.data.frame(region_coherence_clean), 
	function(y) summary(lm(y~ Group+ AgeDeath + Sex + Race, data=pd)))
ptsd_effect_clean = t(sapply(coherence_effects_clean, function(x) x$coef[2,]))
mdd_effect_clean = t(sapply(coherence_effects_clean, function(x) x$coef[3,]))
