###
library(jaffelab)
library(clusterProfiler)
library(dplyr)
library(fields)
library(lattice)
library(RColorBrewer)

dir.create("figures/figure2")

## merged data
load("../../Data/rdas/DEG/AllRegions/geneSet_threeGroups_qSVA_Cortex.rda")
load("../../Data/rdas/DEG/AllRegions/geneSet_threeGroups_qSVA_Amyg.rda")
load("../../Data/rdas/DEG/AllRegions/geneSet_threeGroups_qSVA_allregions.rda")

## by subregion
load("../../Data/rdas/DEG/GOKEGG/dACC/geneSet_threeGroups_qSVA_dACC.rda")
load("../../Data/rdas/DEG/GOKEGG/DLPFC/geneSet_threeGroups_qSVA_DLPFC.rda")
load("../../Data/rdas/DEG/GOKEGG/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg.rda")
load("../../Data/rdas/DEG/GOKEGG/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg.rda")

## combine
goList = list(Cortex = go_cortex,
				Amygdala = go_amyg,
				dACC = go_dacc,
				dlPFC = go_dlpfc,
				MeA = go_medialamyg,
				BLA = go_basoamyg,
				Joint = go_allR)
keggList = list(Cortex = kegg_cortex,
				Amygdala = kegg_amyg,
				dACC = kegg_dacc,
				dlPFC = kegg_dlpfc,
				MeA = kegg_medialamyg,
				BLA = kegg_basoamyg,
				Joint = kegg_allR)
## merge
goList= lapply(goList, as.data.frame)
for(i in seq(along=goList)) {
	goList[[i]]$Region = rep(names(goList)[i])
	goList[[i]]$Cluster = as.character(goList[[i]]$Cluster)
}
keggList= lapply(keggList, as.data.frame)
for(i in seq(along=keggList)) {
	keggList[[i]]$Region = rep(names(keggList)[i])
	keggList[[i]]$Cluster = as.character(keggList[[i]]$Cluster)
	keggList[[i]]$ONTOLOGY = "KEGG"
	keggList[[i]] = keggList[[i]][,colnames(goList[[1]])]
}

## join
go = do.call("rbind", c(goList, keggList))
go$RegionCluster = paste0(go$Region, ":", go$Cluster)

#######################
## filter to mdd #######
#########################

go_mdd = go[go$Cluster %in%
	c("p005all_MDD",	"p005down_MDD",	"p005up_MDD", 
	"BasoAmyg_MDD","BasoAmyg.P.Value_MDD_DOWN","BasoAmyg.P.Value_MDD_UP",
	"MedialAmyg_MDD","MedialAmyg.P.Value_MDD_DOWN","MedialAmyg.P.Value_MDD_UP",
	"DLPFC_MDD","DLPFC.P.Value_MDD_DOWN","DLPFC.P.Value_MDD_UP",
	"dACC_MDD","dACC.P.Value_MDD_DOWN","dACC.P.Value_MDD_UP"),]
go_mdd$Direction = "Both"	
go_mdd$Direction[grep("UP", go_mdd$Cluster,ignore.case=TRUE)] = "Up"	
go_mdd$Direction[grep("DOWN", go_mdd$Cluster,ignore.case=TRUE)] = "Down"	
go_mdd$NewCluster = paste0(go_mdd$Region, "_", go_mdd$Direction)
table(go_mdd$NewCluster )

## convert from long to wide
unique_sets = unique(go_mdd$ID)

go_mdd_wide = do.call("cbind",
	lapply(split(go_mdd, go_mdd$NewCluster), 
		function(x) {
			x[match(unique_sets, x$ID),c("Count", "pvalue", "p.adjust", "geneID")]
		}))
rownames(go_mdd_wide) = unique_sets
go_mdd_wide$ID = unique_sets
go_mdd_wide$Name = go_mdd$Description[match(unique_sets, go_mdd$ID)]
go_mdd_wide$ONTOLOGY = go_mdd$ONTOLOGY[match(unique_sets, go_mdd$ID)]
go_mdd_wide = select(go_mdd_wide, ID,Name, ONTOLOGY, everything())

## and write out table
write.csv(go_mdd_wide, file = "tables/suppTable_SXX_GeneSet_MDD.csv",row.names=FALSE)

#####################
### PTSD vs MDD #####


go_ptsd_vs_mdd = go[go$Cluster %in%
	c("p005all_PTSDvsMDD",	"p005down_PTSDvsMDD",	"p005up_PTSDvsMDD", 
	"BasoAmyg_PTSDvsMDD","BasoAmyg.P.Value_PTSDvsMDD_DOWN","BasoAmyg.P.Value_PTSDvsMDD_UP",
	"MedialAmyg_PTSDvsMDD","MedialAmyg.P.Value_PTSDvsMDD_DOWN","MedialAmyg.P.Value_PTSDvsMDD_UP",
	"DLPFC_PTSDvsMDD","DLPFC.P.Value_PTSDvsMDD_DOWN","DLPFC.P.Value_PTSDvsMDD_UP",
	"dACC_PTSDvsMDD","dACC.P.Value_PTSDvsMDD_DOWN","dACC.P.Value_PTSDvsMDD_UP"),]
go_ptsd_vs_mdd$Direction = "Both"	
go_ptsd_vs_mdd$Direction[grep("UP", go_ptsd_vs_mdd$Cluster,ignore.case=TRUE)] = "Up"	
go_ptsd_vs_mdd$Direction[grep("DOWN", go_ptsd_vs_mdd$Cluster,ignore.case=TRUE)] = "Down"	
go_ptsd_vs_mdd$NewCluster = paste0(go_ptsd_vs_mdd$Region, "_", go_ptsd_vs_mdd$Direction)
table(go_ptsd_vs_mdd$NewCluster )

## convert from long to wide
unique_sets = unique(go_ptsd_vs_mdd$ID)

go_ptsd_vs_mdd_wide = do.call("cbind",
	lapply(split(go_ptsd_vs_mdd, go_ptsd_vs_mdd$NewCluster), 
		function(x) {
			x[match(unique_sets, x$ID),c("Count", "pvalue", "p.adjust", "geneID")]
		}))
rownames(go_ptsd_vs_mdd_wide) = unique_sets
go_ptsd_vs_mdd_wide$ID = unique_sets
go_ptsd_vs_mdd_wide$Name = go_ptsd_vs_mdd$Description[match(unique_sets, go_ptsd_vs_mdd$ID)]
go_ptsd_vs_mdd_wide$ONTOLOGY = go_ptsd_vs_mdd$ONTOLOGY[match(unique_sets, go_ptsd_vs_mdd$ID)]
go_ptsd_vs_mdd_wide = select(go_ptsd_vs_mdd_wide, ID,Name, ONTOLOGY, everything())

## and write out table
write.csv(go_ptsd_vs_mdd_wide, file = "tables/suppTable_SXX_GeneSet_PTSDvsMDD.csv",row.names=FALSE)

############
## plot ####
############

### matrix matrices for heatmaps ###
go_ptsd_vs_mdd_pvalues = go_ptsd_vs_mdd_wide[,grep("pvalue$", colnames(go_ptsd_vs_mdd_wide))]
go_ptsd_vs_mdd_pvalues_up = go_ptsd_vs_mdd_pvalues[,grep("Up", colnames(go_ptsd_vs_mdd_pvalues))]
go_ptsd_vs_mdd_pvalues_down = go_ptsd_vs_mdd_pvalues[,grep("Down", colnames(go_ptsd_vs_mdd_pvalues))]
colnames(go_ptsd_vs_mdd_pvalues_up) = ss(colnames(go_ptsd_vs_mdd_pvalues_up), "_")
colnames(go_ptsd_vs_mdd_pvalues_down) = ss(colnames(go_ptsd_vs_mdd_pvalues_down), "_")

go_ptsd_vs_mdd_fdrs = go_ptsd_vs_mdd_wide[,grep("adjust$", colnames(go_ptsd_vs_mdd_wide))]
go_ptsd_vs_mdd_fdrs_up = go_ptsd_vs_mdd_fdrs[,grep("Up", colnames(go_ptsd_vs_mdd_fdrs))]
go_ptsd_vs_mdd_fdrs_down = go_ptsd_vs_mdd_fdrs[,grep("Down", colnames(go_ptsd_vs_mdd_fdrs))]
colnames(go_ptsd_vs_mdd_fdrs_up) = ss(colnames(go_ptsd_vs_mdd_fdrs_up), "_")
colnames(go_ptsd_vs_mdd_fdrs_down) = ss(colnames(go_ptsd_vs_mdd_fdrs_down), "_")

go_ptsd_vs_mdd_counts = go_ptsd_vs_mdd_wide[,grep("Count", colnames(go_ptsd_vs_mdd_wide))]
go_ptsd_vs_mdd_counts_up = go_ptsd_vs_mdd_counts[,grep("Up", colnames(go_ptsd_vs_mdd_counts))]
go_ptsd_vs_mdd_counts_down = go_ptsd_vs_mdd_counts[,grep("Down", colnames(go_ptsd_vs_mdd_counts))]
colnames(go_ptsd_vs_mdd_counts_up) = ss(colnames(go_ptsd_vs_mdd_counts_up), "_")
colnames(go_ptsd_vs_mdd_counts_down) = ss(colnames(go_ptsd_vs_mdd_counts_down), "_")

#### get some example gene sets
colSums(go_ptsd_vs_mdd_fdrs_up < 0.05, na.rm=TRUE)
colSums(go_ptsd_vs_mdd_fdrs_down < 0.05, na.rm=TRUE)

plot_order = c("Cortex", "dACC", "dlPFC", "Amygdala", "MeA", "BLA")
theSeq = seq(0,8,by=0.01)

########
# DOWN #
lapply(go_ptsd_vs_mdd_fdrs_down[,-6], function(x) {
	xx = go_ptsd_vs_mdd_wide[,c("ID","Name")]
	xx$fdr = x
	xx = xx[which(xx$fdr < 0.05),]
	xx[order(xx$fdr),]
})

down_sets = c("GO:0017080",  "GO:0051592", "GO:0014068", "GO:2001222")
down_pval_plot = go_ptsd_vs_mdd_pvalues_down[down_sets, plot_order]
down_pval_plot[is.na(down_pval_plot)] = 1
rownames(down_pval_plot) = go_ptsd_vs_mdd$Description[match(down_sets,go_ptsd_vs_mdd$ID)]
down_pval_plot[down_pval_plot < 1e-8] = 1e-8
rownames(down_pval_plot)  = gsub("phosphatidylinositol 3-kinase", "PI3K", rownames(down_pval_plot))

pdf("figures/figure4/PTSDvsMDD_geneSet_heatmap_Down.pdf",w=7,h=3.5)
mypal_down = c("white", colorRampPalette(brewer.pal(9,"YlGnBu"))(length(theSeq)))
print(levelplot(-log10(t(down_pval_plot)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=1.6), y=list(cex=1.6)),
		ylab = "", xlab = ""))
dev.off()

########
# up #
lapply(go_ptsd_vs_mdd_fdrs_up[,-6], function(x) {
	xx = go_ptsd_vs_mdd_wide[,c("ID","Name")]
	xx$fdr = x
	xx = xx[which(xx$fdr < 0.05),]
	xx[order(xx$fdr),]
})

up_sets = c("GO:0034702",  "GO:0098982", "GO:0098978", "hsa05033", "GO:0097060")
up_pval_plot = go_ptsd_vs_mdd_pvalues_up[up_sets, plot_order]
up_pval_plot[is.na(up_pval_plot)] = 1
rownames(up_pval_plot) = go_ptsd_vs_mdd$Description[match(up_sets,go_ptsd_vs_mdd$ID)]
up_pval_plot[up_pval_plot < 1e-8] = 1e-8

pdf("figures/figure4/PTSDvsMDD_geneSet_heatmap_Up.pdf",w=5.5,h=3.5)
mypal_up = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(length(theSeq)))
print(levelplot(-log10(t(up_pval_plot)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=1.6), y=list(cex=1.6)),
		ylab = "", xlab = ""))
dev.off()

## and write out table

## look at some for text
go[go$ID == "GO:0048018" & grepl("_PTSD", go$Cluster),]
