
#first, each plot is one brain region with bars colored by GO vs KEGG, disease = PTSD, x is -log10(adj P)
#later, each plot is one brain region with bars colored by disorder (MDD, PTSDvsMDD, onlyPTSD), GO only

####check that i grabbed all incidences of each GO/KEGG term

#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
library(RColorBrewer)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSD.rda")
#load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSD.rda")
#load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD.rda")
#load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD.rda")

#convert to data frame
go_basoamyg_df <- as.data.frame(go_basoamyg)
kegg_basoamyg_df <- as.data.frame(kegg_basoamyg)

#go_basoamyg_ordered <- go_basoamyg_df[order(go_basoamyg_df$p.adjust),]
#kegg_basoamyg_ordered <- kegg_basoamyg_df[order(kegg_basoamyg_df$p.adjust),]

#rename below if not commented out above
go_basoamyg_sig <- go_basoamyg_df[go_basoamyg_df$p.adjust < 0.1,]
kegg_basoamyg_sig <- kegg_basoamyg_df[kegg_basoamyg_df$p.adjust < 0.1,]
go_basoamyg_sig$Type = "GO"
kegg_basoamyg_sig$Type = "KEGG"
kegg_basoamyg_sig$ONTOLOGY = "KEGG"
kegg_basoamyg_sig <- kegg_basoamyg_sig[,c(1,12,2:11)]
both_basoamyg_sig <- rbind(go_basoamyg_sig, kegg_basoamyg_sig)
#still need to incorporate up vs down
#get rid of p005all for this graph
both_basoamyg_sig <- both_basoamyg_sig[both_basoamyg_sig$Cluster == "p005up" | both_basoamyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
both_basoamyg_sig$P.Adj.Calc <- (-log10(both_basoamyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
both_basoamyg_sig$final.P.Adj <- ifelse(both_basoamyg_sig$Cluster == "p005up", both_basoamyg_sig$P.Adj.Calc, -both_basoamyg_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight first pass below)
both_basoamyg_sig2 <- both_basoamyg_sig

pdf("pdf/GOKEGG_BasoAmyg_PTSD.pdf", width=10)
ggplot(data=both_basoamyg_sig2, aes(x = reorder(both_basoamyg_sig2$Description, both_basoamyg_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Type)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-2.5, 2.5, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = both_basoamyg_sig2$Description, y=ifelse(both_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(both_basoamyg_sig2$final.P.Adj < 0, 0, 1), label = both_basoamyg_sig2$Description) +
ggtitle("GO and KEGG BasoAmyg PTSD") 
dev.off()

#change colors?

########################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSD.rda")

#convert to data frame
go_MedialAmyg_df <- as.data.frame(go_medialamyg)
kegg_MedialAmyg_df <- as.data.frame(kegg_medialamyg)

#go_MedialAmyg_ordered <- go_MedialAmyg_df[order(go_MedialAmyg_df$p.adjust),]
#kegg_MedialAmyg_ordered <- kegg_MedialAmyg_df[order(kegg_MedialAmyg_df$p.adjust),]

#rename below if not commented out above
go_MedialAmyg_sig <- go_MedialAmyg_df[go_MedialAmyg_df$p.adjust < 0.1,]
kegg_MedialAmyg_sig <- kegg_MedialAmyg_df[kegg_MedialAmyg_df$p.adjust < 0.1,]
#no kegg sig
go_MedialAmyg_sig$Type = "GO"
#kegg_MedialAmyg_sig$Type = "KEGG"
#kegg_MedialAmyg_sig$ONTOLOGY = "KEGG"
#kegg_MedialAmyg_sig <- kegg_MedialAmyg_sig[,c(1,12,2:11)]
#both_MedialAmyg_sig <- rbind(go_MedialAmyg_sig, kegg_MedialAmyg_sig)
both_MedialAmyg_sig <- go_MedialAmyg_sig
#still need to incorporate up vs down
#get rid of p005all for this graph
both_MedialAmyg_sig <- both_MedialAmyg_sig[both_MedialAmyg_sig$Cluster == "p005up" | both_MedialAmyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
both_MedialAmyg_sig$P.Adj.Calc <- (-log10(both_MedialAmyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
both_MedialAmyg_sig$final.P.Adj <- ifelse(both_MedialAmyg_sig$Cluster == "p005up", both_MedialAmyg_sig$P.Adj.Calc, -both_MedialAmyg_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight first pass below)
both_MedialAmyg_sig2 <- both_MedialAmyg_sig[c(3:6,8,10,11,16,19,20,24,26,29,42,48,53,57,58,61,62,65),]

pdf("pdf/GOKEGG_MedialAmyg_PTSD.pdf", width=15,height=10)
ggplot(data=both_MedialAmyg_sig2, aes(x = reorder(both_MedialAmyg_sig2$Description, both_MedialAmyg_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Type)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-5, 2.5, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = both_MedialAmyg_sig2$Description, y=ifelse(both_MedialAmyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(both_MedialAmyg_sig2$final.P.Adj < 0, 0, 1), label = both_MedialAmyg_sig2$Description) +
ggtitle("GO and KEGG MedialAmyg PTSD") 
dev.off()

#change colors?

########################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD.rda")

#convert to data frame
go_dACC_df <- as.data.frame(go_dACC)
kegg_dACC_df <- as.data.frame(kegg_dACC)

#go_dACC_ordered <- go_dACC_df[order(go_dACC_df$p.adjust),]
#kegg_dACC_ordered <- kegg_dACC_df[order(kegg_dACC_df$p.adjust),]

#rename below if not commented out above
go_dACC_sig <- go_dACC_df[go_dACC_df$p.adjust < 0.1,]
kegg_dACC_sig <- kegg_dACC_df[kegg_dACC_df$p.adjust < 0.1,]
go_dACC_sig$Type = "GO"
kegg_dACC_sig$Type = "KEGG"
kegg_dACC_sig$ONTOLOGY = "KEGG"
kegg_dACC_sig <- kegg_dACC_sig[,c(1,12,2:11)]
both_dACC_sig <- rbind(go_dACC_sig, kegg_dACC_sig)
#still need to incorporate up vs down
#get rid of p005all for this graph
both_dACC_sig <- both_dACC_sig[both_dACC_sig$Cluster == "p005up" | both_dACC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
both_dACC_sig$P.Adj.Calc <- (-log10(both_dACC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
both_dACC_sig$final.P.Adj <- ifelse(both_dACC_sig$Cluster == "p005up", both_dACC_sig$P.Adj.Calc, -both_dACC_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight first pass below)
both_dACC_sig2 <- both_dACC_sig[c(2,3,4,6,7,9,10:14,16,17:21,23),]

pdf("pdf/GOKEGG_dACC_PTSD.pdf", width=10)
ggplot(data=both_dACC_sig2, aes(x = reorder(both_dACC_sig2$Description, both_dACC_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Type)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-2.5, 2.0, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = both_dACC_sig2$Description, y=ifelse(both_dACC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(both_dACC_sig2$final.P.Adj < 0, 0, 1), label = both_dACC_sig2$Description) +
ggtitle("GO and KEGG dACC PTSD") 
dev.off()

#change colors?

########################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD.rda")

#convert to data frame
go_DLPFC_df <- as.data.frame(go_DLPFC)
kegg_DLPFC_df <- as.data.frame(kegg_DLPFC)

#go_DLPFC_ordered <- go_DLPFC_df[order(go_DLPFC_df$p.adjust),]
#kegg_DLPFC_ordered <- kegg_DLPFC_df[order(kegg_DLPFC_df$p.adjust),]

#rename below if not commented out above
go_DLPFC_sig <- go_DLPFC_df[go_DLPFC_df$p.adjust < 0.1,]
kegg_DLPFC_sig <- kegg_DLPFC_df[kegg_DLPFC_df$p.adjust < 0.1,]
go_DLPFC_sig$Type = "GO"
kegg_DLPFC_sig$Type = "KEGG"
kegg_DLPFC_sig$ONTOLOGY = "KEGG"
kegg_DLPFC_sig <- kegg_DLPFC_sig[,c(1,12,2:11)]
both_DLPFC_sig <- rbind(go_DLPFC_sig, kegg_DLPFC_sig)
#still need to incorporate up vs down
#get rid of p005all for this graph
both_DLPFC_sig <- both_DLPFC_sig[both_DLPFC_sig$Cluster == "p005up" | both_DLPFC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
both_DLPFC_sig$P.Adj.Calc <- (-log10(both_DLPFC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
both_DLPFC_sig$final.P.Adj <- ifelse(both_DLPFC_sig$Cluster == "p005up", both_DLPFC_sig$P.Adj.Calc, -both_DLPFC_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight first pass below)
both_DLPFC_sig2 <- both_DLPFC_sig[c(1:5,8:9,14,16:18,23,24,26,27,29,30,33,35,36,37,40,41),]

pdf("pdf/GOKEGG_DLPFC_PTSD.pdf", width=10)
ggplot(data=both_DLPFC_sig2, aes(x = reorder(both_DLPFC_sig2$Description, both_DLPFC_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Type)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-2.5, 2.5, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = both_DLPFC_sig2$Description, y=ifelse(both_DLPFC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(both_DLPFC_sig2$final.P.Adj < 0, 0, 1), label = both_DLPFC_sig2$Description) +
ggtitle("GO and KEGG DLPFC PTSD") 
dev.off()

#change colors?

#add significance line?







