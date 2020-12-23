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

## filter to ptsd
go_ptsd = go[go$Cluster %in%
	c("p005all_PTSD",	"p005down_PTSD",	"p005up_PTSD", 
	"BasoAmyg_PTSD","BasoAmyg.P.Value_PTSD_DOWN","BasoAmyg.P.Value_PTSD_UP",
	"MedialAmyg_PTSD","MedialAmyg.P.Value_PTSD_DOWN","MedialAmyg.P.Value_PTSD_UP",
	"DLPFC_PTSD","DLPFC.P.Value_PTSD_DOWN","DLPFC.P.Value_PTSD_UP",
	"dACC_PTSD","dACC.P.Value_PTSD_DOWN","dACC.P.Value_PTSD_UP"),]
go_ptsd$Direction = "Both"	
go_ptsd$Direction[grep("UP", go_ptsd$Cluster,ignore.case=TRUE)] = "Up"	
go_ptsd$Direction[grep("DOWN", go_ptsd$Cluster,ignore.case=TRUE)] = "Down"	
go_ptsd$NewCluster = paste0(go_ptsd$Region, "_", go_ptsd$Direction)
table(go_ptsd$NewCluster )

## convert from long to wide
unique_sets = unique(go_ptsd$ID)

go_ptsd_wide = do.call("cbind",
	lapply(split(go_ptsd, go_ptsd$NewCluster), 
		function(x) {
			x[match(unique_sets, x$ID),c("Count", "pvalue", "p.adjust", "geneID")]
		}))
rownames(go_ptsd_wide) = unique_sets
go_ptsd_wide$ID = unique_sets
go_ptsd_wide$Name = go_ptsd$Description[match(unique_sets, go_ptsd$ID)]
go_ptsd_wide$ONTOLOGY = go_ptsd$ONTOLOGY[match(unique_sets, go_ptsd$ID)]
go_ptsd_wide = select(go_ptsd_wide, ID,Name, ONTOLOGY, everything())

write.csv(go_ptsd_wide, file = "tables/suppTable_SXX_GeneSet_PTSD.csv",row.names=FALSE)

######################################
##### matrix matrices for heatmaps ###
######################################
go_ptsd_pvalues = go_ptsd_wide[,grep("pvalue$", colnames(go_ptsd_wide))]
go_ptsd_pvalues_up = go_ptsd_pvalues[,grep("Up", colnames(go_ptsd_pvalues))]
go_ptsd_pvalues_down = go_ptsd_pvalues[,grep("Down", colnames(go_ptsd_pvalues))]
colnames(go_ptsd_pvalues_up) = ss(colnames(go_ptsd_pvalues_up), "_")
colnames(go_ptsd_pvalues_down) = ss(colnames(go_ptsd_pvalues_down), "_")

go_ptsd_fdrs = go_ptsd_wide[,grep("adjust$", colnames(go_ptsd_wide))]
go_ptsd_fdrs_up = go_ptsd_fdrs[,grep("Up", colnames(go_ptsd_fdrs))]
go_ptsd_fdrs_down = go_ptsd_fdrs[,grep("Down", colnames(go_ptsd_fdrs))]
colnames(go_ptsd_fdrs_up) = ss(colnames(go_ptsd_fdrs_up), "_")
colnames(go_ptsd_fdrs_down) = ss(colnames(go_ptsd_fdrs_down), "_")

go_ptsd_counts = go_ptsd_wide[,grep("Count", colnames(go_ptsd_wide))]
go_ptsd_counts_up = go_ptsd_counts[,grep("Up", colnames(go_ptsd_counts))]
go_ptsd_counts_down = go_ptsd_counts[,grep("Down", colnames(go_ptsd_counts))]
colnames(go_ptsd_counts_up) = ss(colnames(go_ptsd_counts_up), "_")
colnames(go_ptsd_counts_down) = ss(colnames(go_ptsd_counts_down), "_")

#### get some example gene sets
colSums(go_ptsd_fdrs_up < 0.05, na.rm=TRUE)
colSums(go_ptsd_fdrs_down < 0.05, na.rm=TRUE)

plot_order = c("Cortex", "dACC", "dlPFC", "Amygdala", "MeA", "BLA")
theSeq = seq(0,8,by=0.05)

########
# DOWN #
lapply(go_ptsd_fdrs_down[,-6], function(x) {
	xx = go_ptsd_wide[,c("ID","Name")]
	xx$fdr = x
	xx = xx[which(xx$fdr < 0.01),]
	xx[order(xx$fdr),]
})

down_sets = c("GO:0006959", "GO:0019724", "hsa04610", "GO:0050867",
	"GO:0007159", "GO:0050870", "GO:0048018", 
	"GO:0031012", "GO:0005766", "GO:0098883")
	
down_pval_plot = go_ptsd_pvalues_down[down_sets, plot_order]
down_pval_plot[is.na(down_pval_plot)] = 1
rownames(down_pval_plot) = go_ptsd$Description[match(down_sets,go_ptsd$ID)]
down_pval_plot[down_pval_plot < 1e-8] = 1e-8

pdf("figures/figure2/PTSD_geneSet_heatmap_Down.pdf",w=8.5,h=7)
mypal_down = c("white", colorRampPalette(brewer.pal(9,"YlGnBu"))(length(theSeq)))
print(levelplot(-log10(t(down_pval_plot)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=1.6), y=list(cex=1.6)),
	ylab = "", xlab = "", colorkey=list(labels = list(cex=1.5))))
dev.off()


########
# up #
lapply(go_ptsd_fdrs_up[,-6], function(x) {
	xx = go_ptsd_wide[,c("ID","Name")]
	xx$fdr = x
	xx = xx[which(xx$fdr < 0.05),]
	xx[order(xx$fdr),]
})

up_sets = c("GO:0003682",  "hsa01230", "GO:0016278", 
	"GO:0051931","GO:0017016")
up_pval_plot = go_ptsd_pvalues_up[up_sets, plot_order]
up_pval_plot[is.na(up_pval_plot)] = 1
rownames(up_pval_plot) = go_ptsd$Description[match(up_sets,go_ptsd$ID)]
up_pval_plot[up_pval_plot < 1e-8] = 1e-8

pdf("figures/figure2/PTSD_geneSet_heatmap_Up.pdf",w=8,h=4.5)
mypal_up = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(length(theSeq)))
print(levelplot(-log10(t(up_pval_plot)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=1.6), y=list(cex=1.6)),
		ylab = "", xlab = "",colorkey=list(labels = list(cex=1.5))))
dev.off()
