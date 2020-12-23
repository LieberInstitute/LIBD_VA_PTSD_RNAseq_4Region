
#first, each plot is one brain region with bars colored by GO vs KEGG, disease = PTSD, x is -log10(adj P)
#later, each plot is one brain region with bars colored by disorder (MDD, PTSDvsMDD, onlyPTSD), GO only

####check that i grabbed all incidences of each GO term
#maybe grab ones especially that are opposite between disorders
#BasoAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_onlyPTSD.rda")
go_basoamyg_onlyPTSD <- go_basoamyg
load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_MDD.rda")
go_basoamyg_MDD <- go_basoamyg
load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSDvsMDD.rda")
go_basoamyg_PTSDvsMDD <- go_basoamyg

#change to data frame
go_basoamyg_onlyPTSD_df <- as.data.frame(go_basoamyg_onlyPTSD)
go_basoamyg_MDD_df <- as.data.frame(go_basoamyg_MDD)
go_basoamyg_PTSDvsMDD_df <- as.data.frame(go_basoamyg_PTSDvsMDD)

#keep only significant (adjusted p < 0.1) 
go_basoamyg_onlyPTSD_sig <- go_basoamyg_onlyPTSD_df[go_basoamyg_onlyPTSD_df$p.adjust < 0.1,]
go_basoamyg_MDD_sig <- go_basoamyg_MDD_df[go_basoamyg_MDD_df$p.adjust < 0.1,] #none sig
go_basoamyg_PTSDvsMDD_sig <- go_basoamyg_PTSDvsMDD_df[go_basoamyg_PTSDvsMDD_df$p.adjust < 0.1,]

#uniquely identify the three analyses so can dstinguish after rbinding them
go_basoamyg_onlyPTSD_sig$Analysis = "onlyPTSD"
#go_basoamyg_MDD_sig$Analysis = "MDD" 
go_basoamyg_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

#add in go_basoamyg_MDD_sig if it has sig hits
all_basoamyg_sig <- rbind(go_basoamyg_onlyPTSD_sig, go_basoamyg_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_basoamyg_sig <- all_basoamyg_sig[all_basoamyg_sig$Cluster == "p005up" | all_basoamyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_basoamyg_sig$P.Adj.Calc <- (-log10(all_basoamyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_basoamyg_sig$final.P.Adj <- ifelse(all_basoamyg_sig$Cluster == "p005up", all_basoamyg_sig$P.Adj.Calc, -all_basoamyg_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight
#ex below, still choose final highlighted terms later
all_basoamyg_sig2 <- all_basoamyg_sig[c(2,3,4,6,8,10,12:14,17:19,21,22,31,36,39,44,45:50,52,62,69,72,75,78,80,81),]

pdf("pdf/GO_acrossdisorders_BasoAmyg.pdf", width=10)
ggplot(data=all_basoamyg_sig2, aes(x = reorder(all_basoamyg_sig2$Description, all_basoamyg_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Analysis)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-2.5, 2, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = all_basoamyg_sig2$Description, y=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), label = all_basoamyg_sig2$Description) +
ggtitle("GO BasoAmyg Across Analyses") 
dev.off()

#change colors?

##################################################################

#MedialAmyg
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_onlyPTSD.rda")
go_MedialAmyg_onlyPTSD <- go_medialamyg
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_MDD.rda")
go_MedialAmyg_MDD <- go_medialamyg
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSDvsMDD.rda")
go_MedialAmyg_PTSDvsMDD <- go_medialamyg

#change to data frame
go_MedialAmyg_onlyPTSD_df <- as.data.frame(go_MedialAmyg_onlyPTSD)
go_MedialAmyg_MDD_df <- as.data.frame(go_MedialAmyg_MDD)
go_MedialAmyg_PTSDvsMDD_df <- as.data.frame(go_MedialAmyg_PTSDvsMDD)

#keep only significant (adjusted p < 0.1) 
go_MedialAmyg_onlyPTSD_sig <- go_MedialAmyg_onlyPTSD_df[go_MedialAmyg_onlyPTSD_df$p.adjust < 0.1,]
go_MedialAmyg_MDD_sig <- go_MedialAmyg_MDD_df[go_MedialAmyg_MDD_df$p.adjust < 0.1,]
go_MedialAmyg_PTSDvsMDD_sig <- go_MedialAmyg_PTSDvsMDD_df[go_MedialAmyg_PTSDvsMDD_df$p.adjust < 0.1,]

#uniquely identify the three analyses so can dstinguish after rbinding them
go_MedialAmyg_onlyPTSD_sig$Analysis = "onlyPTSD"
go_MedialAmyg_MDD_sig$Analysis = "MDD" 
go_MedialAmyg_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

all_MedialAmyg_sig <- rbind(go_MedialAmyg_onlyPTSD_sig, go_MedialAmyg_MDD_sig, go_MedialAmyg_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_MedialAmyg_sig <- all_MedialAmyg_sig[all_MedialAmyg_sig$Cluster == "p005up" | all_MedialAmyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_MedialAmyg_sig$P.Adj.Calc <- (-log10(all_MedialAmyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_MedialAmyg_sig$final.P.Adj <- ifelse(all_MedialAmyg_sig$Cluster == "p005up", all_MedialAmyg_sig$P.Adj.Calc, -all_MedialAmyg_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight
#ex below, still choose final highlighted terms later
all_MedialAmyg_sig2 <- all_MedialAmyg_sig[c(1,2,4,5,7,11,14,20,22,23,24,26,28,33,34,38,40,44,47,55,59,62,63,65,77,79,85,86,90,94,103:105),]

pdf("pdf/GO_acrossdisorders_MedialAmyg.pdf", width=10)
ggplot(data=all_MedialAmyg_sig2, aes(x = reorder(all_MedialAmyg_sig2$Description, all_MedialAmyg_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Analysis)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-4.0, 1.5, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab(          "Depleted                      -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = all_MedialAmyg_sig2$Description, y=ifelse(all_MedialAmyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_MedialAmyg_sig2$final.P.Adj < 0, 0, 1), label = all_MedialAmyg_sig2$Description) +
ggtitle("GO MedialAmyg Across Analyses") 
dev.off()

#change colors?

##################################################################

#dACC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_onlyPTSD.rda")
go_dACC_onlyPTSD <- go_dACC
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD.rda")
go_dACC_MDD <- go_dACC
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSDvsMDD.rda")
go_dACC_PTSDvsMDD <- go_dACC

#change to data frame
go_dACC_onlyPTSD_df <- as.data.frame(go_dACC_onlyPTSD)
go_dACC_MDD_df <- as.data.frame(go_dACC_MDD)
go_dACC_PTSDvsMDD_df <- as.data.frame(go_dACC_PTSDvsMDD)

#keep only significant (adjusted p < 0.1) 
go_dACC_onlyPTSD_sig <- go_dACC_onlyPTSD_df[go_dACC_onlyPTSD_df$p.adjust < 0.1,] #none sig
go_dACC_MDD_sig <- go_dACC_MDD_df[go_dACC_MDD_df$p.adjust < 0.1,]
go_dACC_PTSDvsMDD_sig <- go_dACC_PTSDvsMDD_df[go_dACC_PTSDvsMDD_df$p.adjust < 0.1,]



#uniquely identify the three analyses so can dstinguish after rbinding them
go_dACC_onlyPTSD_sig$Analysis = "onlyPTSD"
go_dACC_MDD_sig$Analysis = "MDD" 
go_dACC_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

all_dACC_sig <- rbind(go_dACC_onlyPTSD_sig, go_dACC_MDD_sig, go_dACC_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_dACC_sig <- all_dACC_sig[all_dACC_sig$Cluster == "p005up" | all_dACC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_dACC_sig$P.Adj.Calc <- (-log10(all_dACC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_dACC_sig$final.P.Adj <- ifelse(all_dACC_sig$Cluster == "p005up", all_dACC_sig$P.Adj.Calc, -all_dACC_sig$P.Adj.Calc)



#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight
#ex below, still choose final highlighted terms later
all_dACC_sig2 <- all_dACC_sig[c(1,8,18,24,28,43,46,47,55,64,75,82,143,145,164,182,231,238,266,
	370,384,434,435,484,490,502,504,506,510:512,521:524,530,538,541),]

all_dACC_sig2 <- all_dACC_sig2[c(1:4,5,6,8,11,12,16,17,18,22,23,24,28,33,35,36,37:41),]

pdf("pdf/GO_acrossdisorders_dACC.pdf", width=15,height=10)
ggplot(data=all_dACC_sig2, aes(x = reorder(all_dACC_sig2$Description, all_dACC_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Analysis)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-9.5, 2, by = 1)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = all_dACC_sig2$Description, y=ifelse(all_dACC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_dACC_sig2$final.P.Adj < 0, 0, 1), label = all_dACC_sig2$Description) +
ggtitle("GO dACC Across Analyses") 
dev.off()

#change colors?

##################################################################

#DLPFC
library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_onlyPTSD.rda")
go_DLPFC_onlyPTSD <- go_DLPFC
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD.rda")
go_DLPFC_MDD <- go_DLPFC
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSDvsMDD.rda")
go_DLPFC_PTSDvsMDD <- go_DLPFC

#change to data frame
go_DLPFC_onlyPTSD_df <- as.data.frame(go_DLPFC_onlyPTSD)
go_DLPFC_MDD_df <- as.data.frame(go_DLPFC_MDD)
go_DLPFC_PTSDvsMDD_df <- as.data.frame(go_DLPFC_PTSDvsMDD)

#keep only significant (adjusted p < 0.1) 
go_DLPFC_onlyPTSD_sig <- go_DLPFC_onlyPTSD_df[go_DLPFC_onlyPTSD_df$p.adjust < 0.1,]
go_DLPFC_MDD_sig <- go_DLPFC_MDD_df[go_DLPFC_MDD_df$p.adjust < 0.1,]
go_DLPFC_PTSDvsMDD_sig <- go_DLPFC_PTSDvsMDD_df[go_DLPFC_PTSDvsMDD_df$p.adjust < 0.1,] #none sig

#uniquely identify the three analyses so can dstinguish after rbinding them
go_DLPFC_onlyPTSD_sig$Analysis = "onlyPTSD"
go_DLPFC_MDD_sig$Analysis = "MDD" 
#go_DLPFC_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

#add go_DLPFC_PTSDvsMDD_sig if any sig
all_DLPFC_sig <- rbind(go_DLPFC_onlyPTSD_sig, go_DLPFC_MDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_DLPFC_sig <- all_DLPFC_sig[all_DLPFC_sig$Cluster == "p005up" | all_DLPFC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_DLPFC_sig$P.Adj.Calc <- (-log10(all_DLPFC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_DLPFC_sig$final.P.Adj <- ifelse(all_DLPFC_sig$Cluster == "p005up", all_DLPFC_sig$P.Adj.Calc, -all_DLPFC_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight
#ex below, still choose final highlighted terms later
all_DLPFC_sig2 <- all_DLPFC_sig

pdf("pdf/GO_acrossdisorders_DLPFC.pdf", width=10)
ggplot(data=all_DLPFC_sig2, aes(x = reorder(all_DLPFC_sig2$Description, all_DLPFC_sig2$final.P.Adj), y=as.numeric(final.P.Adj), fill=Analysis)) +
geom_bar(position="dodge",stat="identity") + 
scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
scale_color_brewer(palette="Set2") +
coord_flip() +
ylab( "Depleted            -log10(Adjusted P Value)               Enriched") +
xlab("Description") +
theme(panel.background = element_blank()) +
theme(axis.text.y=element_blank()) +
theme(axis.ticks.y=element_blank()) +
annotate("text", x = all_DLPFC_sig2$Description, y=ifelse(all_DLPFC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_DLPFC_sig2$final.P.Adj < 0, 0, 1), label = all_DLPFC_sig2$Description) +
ggtitle("GO DLPFC Across Analyses") 
dev.off()

#change colors?

#add significance line?
