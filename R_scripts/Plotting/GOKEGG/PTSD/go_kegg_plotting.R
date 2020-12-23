
#MDD

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_MDD.rda")
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_MDD.rda")
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD.rda")
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD.rda")

go_basoamyg_df <- as.data.frame(go_basoamyg)
go_medialamyg_df <- as.data.frame(go_medialamyg)
go_dACC_df <- as.data.frame(go_dACC)
go_DLPFC_df <- as.data.frame(go_DLPFC)
kegg_medialamyg_df <- as.data.frame(kegg_medialamyg)
kegg_basoamyg_df <- as.data.frame(kegg_basoamyg)
kegg_dACC_df <- as.data.frame(kegg_dACC)
kegg_DLPFC_df <- as.data.frame(kegg_DLPFC)
go_basoamyg_ordered <- go_basoamyg_df[order(go_basoamyg_df$p.adjust),]
go_medialamyg_ordered <- go_medialamyg_df[order(go_medialamyg_df$p.adjust),]
go_dACC_ordered <- go_dACC_df[order(go_dACC_df$p.adjust),]
go_DLPFC_ordered <- go_DLPFC_df[order(go_DLPFC_df$p.adjust),]
kegg_basoamyg_ordered <- kegg_basoamyg_df[order(kegg_basoamyg_df$p.adjust),]
kegg_medialamyg_ordered <- kegg_medialamyg_df[order(kegg_medialamyg_df$p.adjust),]
kegg_dACC_ordered <- kegg_dACC_df[order(kegg_dACC_df$p.adjust),]
kegg_DLPFC_ordered <- kegg_DLPFC_df[order(kegg_DLPFC_df$p.adjust),]
go_basoamyg_ordered_topten <- go_basoamyg_ordered[1:10,]
go_medialamyg_ordered_topten <- go_medialamyg_ordered[1:10,]
go_dACC_ordered_topten <- go_dACC_ordered[1:10,]
go_DLPFC_ordered_topten <- go_DLPFC_ordered[1:10,]
kegg_basoamyg_ordered_topten <- kegg_basoamyg_ordered[1:10,]
kegg_medialamyg_ordered_topten <- kegg_medialamyg_ordered[1:10,]
kegg_dACC_ordered_topten <- kegg_dACC_ordered[1:10,]
kegg_DLPFC_ordered_topten <- kegg_DLPFC_ordered[1:10,]


#MDD GO
pdf("pdf/GO_MDD_allregions.pdf", width=8.5)
ggplot(data=go_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO BasoAmyg MDD")

ggplot(data=go_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO MedialAmyg MDD")

ggplot(data=go_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#377eb8") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO dACC MDD")

ggplot(data=go_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
coord_flip() +
ylab("-log10(Adjusted P Value)") +
scale_fill_manual(values=c("#47DA51", "#377eb8")) +
ggtitle("GO DLPFC MDD")

dev.off()

#MDD KEGG 
pdf("pdf/KEGG_MDD_allregions.pdf", width=8.5)
ggplot(data=kegg_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG BasoAmyg MDD")

ggplot(data=kegg_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG MedialAmyg MDD")

ggplot(data=kegg_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG dACC MDD")

ggplot(data=kegg_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG DLPFC MDD")

dev.off()


#PTSD

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSD.rda")
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSD.rda")
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD.rda")
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD.rda")

go_basoamyg_df <- as.data.frame(go_basoamyg)
go_medialamyg_df <- as.data.frame(go_medialamyg)
go_dACC_df <- as.data.frame(go_dACC)
go_DLPFC_df <- as.data.frame(go_DLPFC)
kegg_medialamyg_df <- as.data.frame(kegg_medialamyg)
kegg_basoamyg_df <- as.data.frame(kegg_basoamyg)
kegg_dACC_df <- as.data.frame(kegg_dACC)
kegg_DLPFC_df <- as.data.frame(kegg_DLPFC)
go_basoamyg_ordered <- go_basoamyg_df[order(go_basoamyg_df$p.adjust),]
go_medialamyg_ordered <- go_medialamyg_df[order(go_medialamyg_df$p.adjust),]
go_dACC_ordered <- go_dACC_df[order(go_dACC_df$p.adjust),]
go_DLPFC_ordered <- go_DLPFC_df[order(go_DLPFC_df$p.adjust),]
kegg_basoamyg_ordered <- kegg_basoamyg_df[order(kegg_basoamyg_df$p.adjust),]
kegg_medialamyg_ordered <- kegg_medialamyg_df[order(kegg_medialamyg_df$p.adjust),]
kegg_dACC_ordered <- kegg_dACC_df[order(kegg_dACC_df$p.adjust),]
kegg_DLPFC_ordered <- kegg_DLPFC_df[order(kegg_DLPFC_df$p.adjust),]
go_basoamyg_ordered_topten <- go_basoamyg_ordered[1:10,]
go_medialamyg_ordered_topten <- go_medialamyg_ordered[1:10,]
go_dACC_ordered_topten <- go_dACC_ordered[1:10,]
go_DLPFC_ordered_topten <- go_DLPFC_ordered[1:10,]
kegg_basoamyg_ordered_topten <- kegg_basoamyg_ordered[1:10,]
kegg_medialamyg_ordered_topten <- kegg_medialamyg_ordered[1:10,]
kegg_dACC_ordered_topten <- kegg_dACC_ordered[1:10,]
kegg_DLPFC_ordered_topten <- kegg_DLPFC_ordered[1:10,]


#PTSD GO
pdf("pdf/GO_PTSD_allregions.pdf", width=8.5)
ggplot(data=go_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO BasoAmyg PTSD")

ggplot(data=go_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#377eb8") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO MedialAmyg PTSD")

ggplot(data=go_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO dACC PTSD")

ggplot(data=go_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#377eb8") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO DLPFC PTSD")

dev.off()

#PTSD KEGG 
pdf("pdf/KEGG_PTSD_allregions.pdf", width=8.5)
ggplot(data=kegg_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG BasoAmyg PTSD")

ggplot(data=kegg_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) +
scale_fill_manual(values=c("#47DA51", "#377eb8"))+
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG MedialAmyg PTSD")

ggplot(data=kegg_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#47DA51", "#377eb8"))+
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG dACC PTSD")

ggplot(data=kegg_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG DLPFC PTSD")

dev.off()

#onlyPTSD

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_onlyPTSD.rda")
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_onlyPTSD.rda")
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_onlyPTSD.rda")
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_onlyPTSD.rda")

go_basoamyg_df <- as.data.frame(go_basoamyg)
go_medialamyg_df <- as.data.frame(go_medialamyg)
go_dACC_df <- as.data.frame(go_dACC)
go_DLPFC_df <- as.data.frame(go_DLPFC)
kegg_medialamyg_df <- as.data.frame(kegg_medialamyg)
kegg_basoamyg_df <- as.data.frame(kegg_basoamyg)
kegg_dACC_df <- as.data.frame(kegg_dACC)
kegg_DLPFC_df <- as.data.frame(kegg_DLPFC)
go_basoamyg_ordered <- go_basoamyg_df[order(go_basoamyg_df$p.adjust),]
go_medialamyg_ordered <- go_medialamyg_df[order(go_medialamyg_df$p.adjust),]
go_dACC_ordered <- go_dACC_df[order(go_dACC_df$p.adjust),]
go_DLPFC_ordered <- go_DLPFC_df[order(go_DLPFC_df$p.adjust),]
kegg_basoamyg_ordered <- kegg_basoamyg_df[order(kegg_basoamyg_df$p.adjust),]
kegg_medialamyg_ordered <- kegg_medialamyg_df[order(kegg_medialamyg_df$p.adjust),]
kegg_dACC_ordered <- kegg_dACC_df[order(kegg_dACC_df$p.adjust),]
kegg_DLPFC_ordered <- kegg_DLPFC_df[order(kegg_DLPFC_df$p.adjust),]
go_basoamyg_ordered_topten <- go_basoamyg_ordered[1:10,]
go_medialamyg_ordered_topten <- go_medialamyg_ordered[1:10,]
go_dACC_ordered_topten <- go_dACC_ordered[1:10,]
go_DLPFC_ordered_topten <- go_DLPFC_ordered[1:10,]
kegg_basoamyg_ordered_topten <- kegg_basoamyg_ordered[1:10,]
kegg_medialamyg_ordered_topten <- kegg_medialamyg_ordered[1:10,]
kegg_dACC_ordered_topten <- kegg_dACC_ordered[1:10,]
kegg_DLPFC_ordered_topten <- kegg_DLPFC_ordered[1:10,]


#onlyPTSD GO
pdf("pdf/GO_onlyPTSD_allregions.pdf", width=8.5)
ggplot(data=go_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5,) + 
scale_fill_manual(values=c("#e41a1c", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO BasoAmyg onlyPTSD")

ggplot(data=go_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#e41a1c") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO MedialAmyg onlyPTSD")

ggplot(data=go_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#377eb8") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO dACC onlyPTSD")

ggplot(data=go_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) +
scale_fill_manual(values=c("#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("GO DLPFC onlyPTSD")

dev.off()

#onlyPTSD KEGG 
pdf("pdf/KEGG_onlyPTSD_allregions.pdf", width=8.5)
ggplot(data=kegg_basoamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#e41a1c","#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG BasoAmyg onlyPTSD")

ggplot(data=kegg_medialamyg_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) +
scale_fill_manual(values="#47DA51") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG MedialAmyg onlyPTSD")

ggplot(data=kegg_dACC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values="#377eb8") +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG dACC onlyPTSD")

ggplot(data=kegg_DLPFC_ordered_topten, aes(x = Description, y=-log10(p.adjust), fill=factor(Cluster))) +
geom_bar(position="dodge",stat="identity", width = 0.5) + 
scale_fill_manual(values=c("#47DA51", "#377eb8")) +
coord_flip() +
ylab("-log10(Adjusted P Value)") +
ggtitle("KEGG DLPFC onlyPTSD")

dev.off()

