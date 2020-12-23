
#############################################################################################
##dumbbell charts for GO across disorders
#Style info
#note, only including PTSD analysis for these figures if overlapping
#############################################################################################

#PTSDvsMDD: #E3712B
#onlyPTSD: #3881B0
#MDD: #CB8CAD
#PTSD: #000000
#cool green: #42892E
#other green: #56A255
#other green: #52A55E 
#light green: #85C08E
#redish color: #D13939
#light red: #D67272

#color selection
library(RColorBrewer)
brewer.pal(9, "Set1")
display.brewer.pal(9, "Set1")

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getPalette(16)
# "#E41A1C" "#874F6F" "#3881B0" "#449B75" "#56A255" "#7E6E85" "#AC5782"
# "#E3712B" "#FFA10D" "#FFE528" "#E1C62F" "#B16C29" "#C66764" "#F17EB4"
# "#CB8CAD" "#999999"

plot(rep(1,16),col=getPalette(16),pch=19,cex=3)

#############################################################################################
#BasoAmyg
#############################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(RColorBrewer)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSD.rda")
go_basoamyg_PTSD <- go_basoamyg
load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_onlyPTSD.rda")
go_basoamyg_onlyPTSD <- go_basoamyg
load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_MDD.rda")
go_basoamyg_MDD <- go_basoamyg
load("rdas/BasoAmyg/geneSet_threeGroups_qSVA_BasoAmyg_PTSDvsMDD.rda")
go_basoamyg_PTSDvsMDD <- go_basoamyg

#change to data frame
go_basoamyg_PTSD_df <- as.data.frame(go_basoamyg_PTSD)
go_basoamyg_onlyPTSD_df <- as.data.frame(go_basoamyg_onlyPTSD)
go_basoamyg_MDD_df <- as.data.frame(go_basoamyg_MDD)
go_basoamyg_PTSDvsMDD_df <- as.data.frame(go_basoamyg_PTSDvsMDD)

#keep only significant (adjusted p < 0.1)
go_basoamyg_PTSD_sig <- go_basoamyg_PTSD_df[go_basoamyg_PTSD_df$p.adjust < 0.1,] 
go_basoamyg_onlyPTSD_sig <- go_basoamyg_onlyPTSD_df[go_basoamyg_onlyPTSD_df$p.adjust < 0.1,]
go_basoamyg_MDD_sig <- go_basoamyg_MDD_df[go_basoamyg_MDD_df$p.adjust < 0.1,] #none sig
go_basoamyg_PTSDvsMDD_sig <- go_basoamyg_PTSDvsMDD_df[go_basoamyg_PTSDvsMDD_df$p.adjust < 0.1,]

#uniquely identify the three analyses so can dstinguish after rbinding them
go_basoamyg_onlyPTSD_sig$Analysis = "onlyPTSD"
#go_basoamyg_MDD_sig$Analysis = "MDD" 
go_basoamyg_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"
go_basoamyg_PTSD_sig$Analysis = "PTSD"

###here
#add in go_basoamyg_MDD_sig if it has sig hits
all_basoamyg_sig <- rbind(go_basoamyg_PTSD_sig,go_basoamyg_onlyPTSD_sig, go_basoamyg_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_basoamyg_sig <- all_basoamyg_sig[all_basoamyg_sig$Cluster == "p005up" | all_basoamyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_basoamyg_sig$P.Adj.Calc <- (-log10(all_basoamyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_basoamyg_sig$final.P.Adj <- ifelse(all_basoamyg_sig$Cluster == "p005up", all_basoamyg_sig$P.Adj.Calc, -all_basoamyg_sig$P.Adj.Calc)


#note this includes all three kinds of GO terms (BP, MF, CC)
#select descriptions we want to highlight
#prioritize descriptions in common as well as higher gene ratio


all_basoamyg_sig[which(all_basoamyg_sig$Description %in% 
	intersect(all_basoamyg_sig$Description[all_basoamyg_sig$Analysis == "PTSD"],intersect(all_basoamyg_sig$Description[all_basoamyg_sig$Analysis=="onlyPTSD"],
	all_basoamyg_sig$Description[all_basoamyg_sig$Analysis=="PTSDvsMDD"]))),]

intersect(all_basoamyg_sig$Description[all_basoamyg_sig$Analysis == "PTSDvsMDD"],all_basoamyg_sig$Description[all_basoamyg_sig$Analysis=="onlyPTSD"])


#turn gene ratio into decimal so we can order on it
all_basoamyg_sig$GeneRatioDec <- sapply(all_basoamyg_sig$GeneRatio, function(x) eval(parse(text = x)))
all_basoamyg_sig <- all_basoamyg_sig[order(all_basoamyg_sig$GeneRatioDec, decreasing=TRUE),]

all_basoamyg_sig[all_basoamyg_sig$Analysis == "onlyPTSD",][1:20,]
all_basoamyg_sig[all_basoamyg_sig$Analysis == "PTSDvsMDD",][1:20,]
all_basoamyg_sig[all_basoamyg_sig$Analysis == "PTSD",]


basoamyg_keep <- c("regulation of glutamate receptor signaling pathway","regulation of neurotransmitter receptor activity","learning","positive regulation of phosphatidylinositol 3-kinase signaling",
	"regulation of signaling receptor activity","behavior","presynapse","cognition","regulation of membrane potential","action potential","G protein-coupled receptor activity","inhibitory synapse",
	"synaptic transmission, GABAergic","negative regulation of macrophage activation","positive regulation of protein localization to cell surface","response to pain",
	"neuron migration", "negative regulation of neuroinflammatory response","glutamate binding","mismatched DNA binding","neuropeptide signaling pathway","long-term memory",
	"multicellular organismal response to stress","regulation of NMDA receptor activity","spliceosomal snRNP complex","lysophospholipase activity")

all_basoamyg_sig2 <- all_basoamyg_sig[which(all_basoamyg_sig$Description %in% basoamyg_keep),]

#dumbbell

library("ggalt")
library("tidyr")
#library("bbplot")
#library(gridExtra)
library(dplyr)
library(stringr)

font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("onlyPTSD" = "#3881B0","PTSDvsMDD" = "#E3712B")

all_basoamyg_sig2<- all_basoamyg_sig2[order(all_basoamyg_sig2$final.P.Adj,decreasing=FALSE),]

all_basoamyg_sig2$Description <- factor(all_basoamyg_sig2$Description, levels = unique(all_basoamyg_sig2$Description))

##weirdness so that ggplot2 layers work out correctly: in the plots below, I added a specific dumbbell layer to get the semgments to show correctly

##for pdf plots made on the cluster, i need position_nudge just due to graphics display compatibility issues

##The below generates a plot with the point/dumbbell color set by direction logFC and the line/segment color set by analysis)
#can align points with segment
pdf("pdf/GO_BasoAmyg_AcrossDisorders_segmentanddumbbell_differentcolors.pdf",width=10)
ggplot(data=all_basoamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description, colour=Analysis)) +
geom_dumbbell(colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=3) +
geom_dumbbell(data=subset(all_basoamyg_sig2,rownames(all_basoamyg_sig2) == "2567"), colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=3) +
geom_point(data=subset(all_basoamyg_sig2,all_basoamyg_sig2$final.P.Adj > 0),aes(y = Description, x=final.P.Adj),color="#52A55E",size=4.25, position=position_nudge(y=0.05)) + 
geom_point(data=subset(all_basoamyg_sig2,all_basoamyg_sig2$final.P.Adj < 0),aes(y = Description, x=final.P.Adj),color="#D13939",size=4.25,position=position_nudge(y=0.05)) +
scale_colour_manual(values=analysis_colors) +
guides(colour = guide_legend(override.aes = list(shape = 15))) +
annotate("text", y = all_basoamyg_sig2$Description, x=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_basoamyg_sig2$Description,50),vjust=0.5, lineheight = 0.85) +
ggtitle("GO BasoAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2.5,2.5),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)) , breaks=seq(-2.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=20,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=16,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=16,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=14,color=text_color),
	axis.title = element_text(family=font, size=16,color=text_color,face="bold", margin=ggplot2::margin(9,0,1,0)),
	axis.text = element_text(family=font, size=14,color=text_color),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()


##The below generates a plot with the dumbbell end colored by Analysis and the segment is a solid grey line
#winner winner chicken dinner

#get x axis tick colors, code depends on axis limits
x_colors<- c(rep("#C24A4A",5),"#222222",rep("#5A9C64",5))

pdf("pdf/GO_BasoAmyg_AcrossDisorders_dumbbell_latest.pdf",width=10)
ggplot(data=all_basoamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_dumbbell(data=subset(all_basoamyg_sig2,rownames(all_basoamyg_sig2) == "2567"), colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) +
geom_point(data=all_basoamyg_sig2,aes(y=all_basoamyg_sig2$Description, x=all_basoamyg_sig2$final.P.Adj,color=Analysis),size=4.25, position=position_nudge(y=0.05)) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_basoamyg_sig2$Description, x=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_basoamyg_sig2$Description,50),vjust=0.5, lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO BasoAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2.5,2.5),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-2.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()


##The below generates a plot with the segment and dumbbell end colored by Analysis
#i added a border to the dumbbell for fun/visual interest but remove border from legend
pdf("pdf/GO_BasoAmyg_AcrossDisorders_segmentanddumbbell_samecolors.pdf",width=10)
ggplot(data=all_basoamyg_sig2, aes(x=0, xend=final.P.Adj,y =Description, colour=Analysis, fill=Analysis)) +
geom_dumbbell(colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=4) + 
geom_dumbbell(data=subset(all_basoamyg_sig2,rownames(all_basoamyg_sig2) == "2567"), colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=4) +
geom_point(data=all_basoamyg_sig2,aes(y=all_basoamyg_sig2$Description, x=all_basoamyg_sig2$final.P.Adj),size=4, position=position_nudge(y=0.05),shape=16) +
geom_point(data=all_basoamyg_sig2,aes(y=all_basoamyg_sig2$Description, x=all_basoamyg_sig2$final.P.Adj),size=4, position=position_nudge(y=0.05), colour="#222222", shape=21, stroke=0.75, show.legend=FALSE) +
scale_colour_manual(values=analysis_colors) +
scale_fill_manual(values=analysis_colors) +
annotate("text", y = all_basoamyg_sig2$Description, x=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_basoamyg_sig2$Description,50),vjust=0.5, lineheight = 0.85) +
ggtitle("GO BasoAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2.5,2.5),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-2.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=20,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=16,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=16,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=14,color=text_color),
	axis.title = element_text(family=font, size=16,color=text_color,face="bold", margin=ggplot2::margin(9,0,1,0)),
	axis.text = element_text(family=font, size=14,color=text_color),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()


#png
##The below generates a plot with the dumbbell end colored by Analysis and the segment is a solid grey line
#i have to generate this on my laptop, not the cluster, since ggsave has font issues with the cluster

##make size of dot relative to generatio (keep consistent for all 4 regions)
#bin values so starting value is min and you fall into that bin until you hit the next cutoff pt, no rounding
#0.01, 0.05, 0.1, 0.15, 0.2 --> point sizes: 4, 5.5, 7, 8.5

all_basoamyg_sig2$point_sizes<- ifelse(all_basoamyg_sig2$GeneRatioDec < 0.01, NA,ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.01 & all_basoamyg_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.05 & all_basoamyg_sig2$GeneRatioDec < 0.1, 6, ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.1 & all_basoamyg_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_basoamyg_sig2$GeneRatioDec >= 0.15 & all_basoamyg_sig2$GeneRatioDec < 0.2, 10, NA)))))
#save object so that i can make figure as png on my laptop (off the cluster)
save(all_basoamyg_sig2, file="rdas/GO_BasoAmyg_p.adjust0.1_plotinput.rda")



#eventually add symbols to delineate ontologies...
#circle = CC
#square = BP
#triangle = MF
#all_basoamyg_sig2$Symbols <- ifelse(all_basoamyg_sig2$ONTOLOGY == "CC",1, ifelse(all_basoamyg_sig2$ONTOLOGY == "BP",0,2))  

x_colors<- c(rep("#C24A4A",5),"#222222",rep("#5A9C64",5))
font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("onlyPTSD" = "#3881B0","PTSDvsMDD" = "#E3712B")
point_sizes<-all_basoamyg_sig2$point_sizes 										

										
basoamyg_dumbell <- ggplot(data=all_basoamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_dumbbell(data=subset(all_basoamyg_sig2,rownames(all_basoamyg_sig2) == "2567"), colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) +
geom_point(data=all_basoamyg_sig2,aes(y=all_basoamyg_sig2$Description, x=all_basoamyg_sig2$final.P.Adj,color=Analysis),size=point_sizes) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_basoamyg_sig2$Description, x=ifelse(all_basoamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_basoamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_basoamyg_sig2$Description,70),color=text_color, family=font, vjust=0.5,lineheight = 0.85) +
ggtitle("GO BasoAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2.5,2.5),labels = c(seq(2.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-2.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=basoamyg_dumbell,file="png/GO_BasoAmyg_AcrossDisorders_dumbbell_latest.png",width=10,height=8.5,dpi=300)


#############################################################################################
#MedialAmyg
#############################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(RColorBrewer)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSD.rda",verbose=TRUE)
go_MedialAmyg_PTSD <- go_medialamyg
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_onlyPTSD.rda")
go_MedialAmyg_onlyPTSD <- go_medialamyg
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_MDD.rda")
go_MedialAmyg_MDD <- go_medialamyg
load("rdas/MedialAmyg/geneSet_threeGroups_qSVA_MedialAmyg_PTSDvsMDD.rda")
go_MedialAmyg_PTSDvsMDD <- go_medialamyg

#change to data frame
go_MedialAmyg_PTSD_df <- as.data.frame(go_MedialAmyg_PTSD)
go_MedialAmyg_onlyPTSD_df <- as.data.frame(go_MedialAmyg_onlyPTSD)
go_MedialAmyg_MDD_df <- as.data.frame(go_MedialAmyg_MDD)
go_MedialAmyg_PTSDvsMDD_df <- as.data.frame(go_MedialAmyg_PTSDvsMDD)

#keep only significant (adjusted p < 0.1)
go_MedialAmyg_PTSD_sig <- go_MedialAmyg_PTSD_df[go_MedialAmyg_PTSD_df$p.adjust < 0.1,] 
go_MedialAmyg_onlyPTSD_sig <- go_MedialAmyg_onlyPTSD_df[go_MedialAmyg_onlyPTSD_df$p.adjust < 0.1,]
go_MedialAmyg_MDD_sig <- go_MedialAmyg_MDD_df[go_MedialAmyg_MDD_df$p.adjust < 0.1,]
go_MedialAmyg_PTSDvsMDD_sig <- go_MedialAmyg_PTSDvsMDD_df[go_MedialAmyg_PTSDvsMDD_df$p.adjust < 0.1,]

#uniquely identify the three analyses so can dstinguish after rbinding them
go_MedialAmyg_PTSD_sig$Analysis = "PTSD"
go_MedialAmyg_onlyPTSD_sig$Analysis = "onlyPTSD"
go_MedialAmyg_MDD_sig$Analysis = "MDD" 
go_MedialAmyg_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

all_MedialAmyg_sig <- rbind(go_MedialAmyg_PTSD_sig, go_MedialAmyg_onlyPTSD_sig, go_MedialAmyg_MDD_sig, go_MedialAmyg_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_MedialAmyg_sig <- all_MedialAmyg_sig[all_MedialAmyg_sig$Cluster == "p005up" | all_MedialAmyg_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_MedialAmyg_sig$P.Adj.Calc <- (-log10(all_MedialAmyg_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_MedialAmyg_sig$final.P.Adj <- ifelse(all_MedialAmyg_sig$Cluster == "p005up", all_MedialAmyg_sig$P.Adj.Calc, -all_MedialAmyg_sig$P.Adj.Calc)

all_medialamyg_sig <- all_MedialAmyg_sig

all_medialamyg_sig[which(all_medialamyg_sig$Description %in% 
	intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis == "PTSD"],intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="onlyPTSD"],
	intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="MDD"],
	all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="PTSDvsMDD"])))),]

intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis == "PTSDvsMDD"],all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="onlyPTSD"])

all_medialamyg_sig$Description[which(all_medialamyg_sig$Description %in% 
	intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis == "PTSDvsMDD"],intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis == "PTSD"],
	all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="onlyPTSD"])))]

intersect(all_medialamyg_sig$Description[all_medialamyg_sig$Analysis == "MDD"],all_medialamyg_sig$Description[all_medialamyg_sig$Analysis=="PTSD"])
                                           

all_medialamyg_sig$GeneRatioDec <- sapply(all_medialamyg_sig$GeneRatio, function(x) eval(parse(text = x)))

all_medialamyg_sig <- all_medialamyg_sig[order(all_medialamyg_sig$GeneRatioDec, decreasing=TRUE),]

all_medialamyg_sig[all_medialamyg_sig$Analysis == "onlyPTSD",][1:20,"Description"]
all_medialamyg_sig[all_medialamyg_sig$Analysis == "PTSDvsMDD",][1:20,"Description"]
all_medialamyg_sig[all_medialamyg_sig$Analysis == "PTSD",][1:20,"Description"]
all_medialamyg_sig[all_medialamyg_sig$Analysis == "MDD",][1:20,"Description"]

medialamyg_keep <- c("cation transmembrane transporter activity","extracellular matrix organization","regulation of signaling receptor activity",
	"glucose metabolic process", "ATP-dependent helicase activity","passive transmembrane transporter activity","proximal promoter sequence-specific DNA binding","voltage-gated ion channel activity",
	"regulation of behavior","ligand-gated ion channel activity", "regulation of signaling receptor activity","mitochondrial transport","cellular lipid catabolic process",
	"cholesterol metabolic process", "collagen fibril organization", "fatty acid transmembrane transport", "dopaminergic neuron differentiation", "platelet-derived growth factor binding",
	"transmembrane receptor protein serine/threonine kinase signaling pathway","histone deacetylase binding","telomeric DNA binding","autophagy","cellular ketone metabolic process")

all_medialamyg_sig2 <- all_medialamyg_sig[which(all_medialamyg_sig$Description %in% medialamyg_keep),]

#dumbbell

library("ggalt")
library("tidyr")
library(dplyr)
library(stringr)

font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("MDD" = "#CB8CAD", "PTSD" = "#000000","onlyPTSD" = "#3881B0","PTSDvsMDD" = "#E3712B")
x_colors<- c(rep("#C24A4A",11),"#222222",rep("#5A9C64",4))

all_medialamyg_sig2<- all_medialamyg_sig2[order(all_medialamyg_sig2$final.P.Adj,decreasing=FALSE),]

all_medialamyg_sig2$Description <- factor(all_medialamyg_sig2$Description, levels = unique(all_medialamyg_sig2$Description))
all_medialamyg_sig2$Analysis <- factor(all_medialamyg_sig2$Analysis, levels = c("MDD","PTSD","onlyPTSD","PTSDvsMDD"))

all_medialamyg_sig2$point_sizes<- ifelse(all_medialamyg_sig2$GeneRatioDec < 0.01, NA,ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.01 & all_medialamyg_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.05 & all_medialamyg_sig2$GeneRatioDec < 0.1, 6, ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.1 & all_medialamyg_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.15 & all_medialamyg_sig2$GeneRatioDec < 0.2, 10, NA)))))

#save object so that i can make figure as png on my laptop (off the cluster)
save(all_medialamyg_sig2, file="rdas/GO_MedialAmyg_p.adjust0.1_plotinput.rda")

##The below generates a plot with the dumbbell end colored by Analysis and the segment is a solid grey line

wrap_MeA <- stringr::str_wrap(all_medialamyg_sig2$Description,55)
wrap_MeA[7] <- "transmembrane receptor protein serine/\nthreonine kinase signaling pathway"

pdf("pdf/GO_MedialAmyg_AcrossDisorders_dumbbell_latest.pdf",width=10)
ggplot(data=all_medialamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_medialamyg_sig2,aes(y=all_medialamyg_sig2$Description, x=all_medialamyg_sig2$final.P.Adj,color=Analysis),size=point_sizes, position=position_nudge(y=0.05)) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_medialamyg_sig2$Description, x=ifelse(all_medialamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_medialamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_medialamyg_sig2$Description,45),vjust=0.5, lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO MedialAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-5.5,2),labels = c(seq(5.5, 0, by = -0.5),seq(0.5, 2, by = 0.5)), breaks=seq(-5.5, 2, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin = margin(r=-20)),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()



#png
#i have to generate this on my laptop, not the cluster, since ggsave has font issues with the cluster
pointsizes <- c("0.05" = 4, "0.1" = 6,"0.15" = 8,"0.2" = 10)

all_medialamyg_sig2$point_sizes<- ifelse(all_medialamyg_sig2$GeneRatioDec < 0.01, NA,ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.01 & all_medialamyg_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.05 & all_medialamyg_sig2$GeneRatioDec < 0.1, 6, ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.1 & all_medialamyg_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_medialamyg_sig2$GeneRatioDec >= 0.15 & all_medialamyg_sig2$GeneRatioDec < 0.2, 10, NA)))))

all_medialamyg_sig2$Size <- ifelse(all_medialamyg_sig2$point_sizes == 4, 0.05, ifelse(all_medialamyg_sig2$point_sizes == 6, 0.1, 
	ifelse(all_medialamyg_sig2$point_sizes == 8, 0.15, ifelse(all_medialamyg_sig2$point_sizes == 10, 0.2, NA))))

all_medialamyg_sig2$Size <-factor(all_medialamyg_sig2$Size, levels = c("0.05","0.1","0.15","0.2"))


wrap_MeA <- stringr::str_wrap(all_medialamyg_sig2$Description,55)
wrap_MeA[7] <- "transmembrane receptor protein serine/\nthreonine kinase signaling pathway"

medialamyg_dumbell <- ggplot(data=all_medialamyg_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_medialamyg_sig2,aes(y=all_medialamyg_sig2$Description, x=all_medialamyg_sig2$final.P.Adj,color=Analysis,size=Size)) +
scale_colour_manual(values=analysis_colors) +
scale_size_manual(values=pointsizes) +
annotate("text", y = all_medialamyg_sig2$Description, x=ifelse(all_medialamyg_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_medialamyg_sig2$final.P.Adj < 0, 0, 1), 
	label = wrap_MeA,vjust=0.5, lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO MedialAmyg", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-5.5,2),labels = c(seq(5.5, 0, by = -0.5),seq(0.5, 2, by = 0.5)), breaks=seq(-5.5, 2, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.77),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(t=0.25,r=0,b=9,l=0),hjust = 0.77),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5), hjust=0.79),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin = margin(r=-20)),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=medialamyg_dumbell,file="png/GO_MedialAmyg_AcrossDisorders_dumbbell_test.png",width=10,height=8.5,dpi=300)


#############################################################################################
#dACC
#############################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(RColorBrewer)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSD.rda")
go_dACC_PTSD <- go_dACC
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_onlyPTSD.rda")
go_dACC_onlyPTSD <- go_dACC
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_MDD.rda")
go_dACC_MDD <- go_dACC
load("rdas/dACC/geneSet_threeGroups_qSVA_dACC_PTSDvsMDD.rda")
go_dACC_PTSDvsMDD <- go_dACC

#change to data frame
go_dACC_PTSD_df <- as.data.frame(go_dACC_PTSD)
go_dACC_onlyPTSD_df <- as.data.frame(go_dACC_onlyPTSD)
go_dACC_MDD_df <- as.data.frame(go_dACC_MDD)
go_dACC_PTSDvsMDD_df <- as.data.frame(go_dACC_PTSDvsMDD)

#keep only significant (adjusted p < 0.1)
go_dACC_PTSD_sig <- go_dACC_PTSD_df[go_dACC_PTSD_df$p.adjust < 0.1,] 
go_dACC_onlyPTSD_sig <- go_dACC_onlyPTSD_df[go_dACC_onlyPTSD_df$p.adjust < 0.1,] #none sig
go_dACC_MDD_sig <- go_dACC_MDD_df[go_dACC_MDD_df$p.adjust < 0.1,]
go_dACC_PTSDvsMDD_sig <- go_dACC_PTSDvsMDD_df[go_dACC_PTSDvsMDD_df$p.adjust < 0.1,]

#uniquely identify the three analyses so can dstinguish after rbinding them
go_dACC_PTSD_sig$Analysis = "PTSD"
#go_dACC_onlyPTSD_sig$Analysis = "onlyPTSD"
go_dACC_MDD_sig$Analysis = "MDD" 
go_dACC_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

all_dACC_sig <- rbind(go_dACC_PTSD_sig, go_dACC_MDD_sig, go_dACC_PTSDvsMDD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_dACC_sig <- all_dACC_sig[all_dACC_sig$Cluster == "p005up" | all_dACC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_dACC_sig$P.Adj.Calc <- (-log10(all_dACC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_dACC_sig$final.P.Adj <- ifelse(all_dACC_sig$Cluster == "p005up", all_dACC_sig$P.Adj.Calc, -all_dACC_sig$P.Adj.Calc)

#all_dACC_sig[which(all_dACC_sig$Description %in% 
#	intersect(all_dACC_sig$Description[all_dACC_sig$Analysis == "PTSD"],intersect(all_dACC_sig$Description[all_dACC_sig$Analysis=="onlyPTSD"],
#	intersect(all_dACC_sig$Description[all_dACC_sig$Analysis=="MDD"],
#	all_dACC_sig$Description[all_dACC_sig$Analysis=="PTSDvsMDD"])))),]

intersect(all_dACC_sig$Description[all_dACC_sig$Analysis == "PTSDvsMDD"],all_dACC_sig$Description[all_dACC_sig$Analysis=="PTSD"])

all_dACC_sig[which(all_dACC_sig$Description %in% 
	intersect(all_dACC_sig$Description[all_dACC_sig$Analysis == "PTSDvsMDD"],intersect(all_dACC_sig$Description[all_dACC_sig$Analysis == "PTSD"],
	all_dACC_sig$Description[all_dACC_sig$Analysis=="MDD"]))),]

intersect(all_dACC_sig$Description[all_dACC_sig$Analysis == "MDD"],all_dACC_sig$Description[all_dACC_sig$Analysis=="PTSD"])                                         

all_dACC_sig$GeneRatioDec <- sapply(all_dACC_sig$GeneRatio, function(x) eval(parse(text = x)))

all_dACC_sig <- all_dACC_sig[order(all_dACC_sig$GeneRatioDec, decreasing=TRUE),]


#all_dACC_sig[all_dACC_sig$Analysis == "onlyPTSD",][1:20,"Description"]
all_dACC_sig[all_dACC_sig$Analysis == "PTSDvsMDD",][1:20,"Description"]
all_dACC_sig[all_dACC_sig$Analysis == "PTSD",][1:20,"Description"]
all_dACC_sig[all_dACC_sig$Analysis == "MDD",][1:20,"Description"]

dACC_keep <- c("plasma membrane protein complex","SH2 domain binding","humoral immune response","complement activation","acute inflammatory response",
	"positive regulation of smooth muscle contraction","transport vesicle","postsynaptic membrane","asymmetric synapse","neuron to neuron synapse",
	"transmembrane transporter complex","distal axon","cation transmembrane transporter activity","microglial cell activation","neuroinflammatory response",
	"toll-like receptor signaling pathway","negative regulation of cell-cell adhesion","voltage-gated ion channel activity","exocytic vesicle membrane","glutamatergic synapse",
	"cholesterol binding","insulin-like growth factor receptor binding","lymphocyte activation","synapse pruning","adaptive immune response","regulation of vesicle-mediated transport")

all_dACC_sig2 <- all_dACC_sig[which(all_dACC_sig$Description %in% dACC_keep),]

#dumbbell

library("ggalt")
library("tidyr")
library(dplyr)
library(stringr)

font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("MDD" = "#CB8CAD", "PTSD" = "#000000","PTSDvsMDD" = "#E3712B")
x_colors<- c(rep("#C24A4A",10),"#222222",rep("#5A9C64",5))

all_dACC_sig2<- all_dACC_sig2[order(all_dACC_sig2$final.P.Adj,decreasing=FALSE),]

all_dACC_sig2$Description <- factor(all_dACC_sig2$Description, levels = unique(all_dACC_sig2$Description))

all_dACC_sig2$point_sizes<- ifelse(all_dACC_sig2$GeneRatioDec < 0.01, NA,ifelse(all_dACC_sig2$GeneRatioDec >= 0.01 & all_dACC_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_dACC_sig2$GeneRatioDec >= 0.05 & all_dACC_sig2$GeneRatioDec < 0.1, 6, ifelse(all_dACC_sig2$GeneRatioDec >= 0.1 & all_dACC_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_dACC_sig2$GeneRatioDec >= 0.15 & all_dACC_sig2$GeneRatioDec < 0.2, 10, NA)))))

#save object so that i can make figure as png on my laptop (off the cluster)
save(all_dACC_sig2, file="rdas/GO_dACC_p.adjust0.1_plotinput.rda")

##The below generates a plot with the dumbbell end colored by Analysis and the segment is a solid grey line

font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("MDD" = "#CB8CAD", "PTSD" = "#000000","PTSDvsMDD" = "#E3712B")
x_colors<- c(rep("#C24A4A",21),"#222222",rep("#5A9C64",5))

justify_x <- ifelse(all_dACC_sig2$Description=="plasma membrane protein complex", 0, ifelse(all_dACC_sig2$Description=="SH2 domain binding", 0,
	ifelse(all_dACC_sig2$final.P.Adj < 0,0.1,-0.1)))

justify_h <- ifelse(all_dACC_sig2$Description=="plasma membrane protein complex" | all_dACC_sig2$Description=="SH2 domain binding", 0.5,
	ifelse(all_dACC_sig2$final.P.Adj < 0, 0, 1))

wrap_dACC<- stringr::str_wrap(all_dACC_sig2$Description,45)
wrap_dACC[c(14,18)] <- "positive regulation of smooth \nmuscle contraction"

pdf("pdf/GO_dACC_AcrossDisorders_dumbbell_latest.pdf",width=10)
ggplot(data=all_dACC_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_dACC_sig2,aes(y=all_dACC_sig2$Description, x=all_dACC_sig2$final.P.Adj,color=Analysis),size=point_sizes, position=position_nudge(y=0.05)) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_dACC_sig2$Description, , x=justify_x, hjust = justify_h, 
	label = wrap_dACC, lineheight = 0.85,color=text_color, family=font, vjust=0.5) +
ggtitle("GO dACC", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-10.5,2.5),labels = c(seq(10.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-10.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin = margin(r=-50)),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()

#png
#i have to generate this on my laptop, not the cluster, since ggsave has font issues with the cluster
load("GO_dACC_p.adjust0.1_plotinput.rda")
font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("MDD" = "#CB8CAD", "PTSD" = "#000000","PTSDvsMDD" = "#E3712B")
x_colors<- c(rep("#C24A4A",21),"#222222",rep("#5A9C64",5))
point_sizes<-all_dACC_sig2$point_sizes 										

justify_x <- ifelse(all_dACC_sig2$final.P.Adj < 0,0.1,-0.1)

justify_h <- ifelse(all_dACC_sig2$final.P.Adj < 0, 0, 1)

wrap_dACC<- stringr::str_wrap(all_dACC_sig2$Description,45)
wrap_dACC[c(14,18)] <- "positive regulation of smooth \nmuscle contraction"


dACC_dumbell <- ggplot(data=all_dACC_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_dACC_sig2,aes(y=all_dACC_sig2$Description, x=all_dACC_sig2$final.P.Adj,color=Analysis),size=point_sizes) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_dACC_sig2$Description[c(-23,-10,-8,-32)], x=justify_x[c(-23,-10,-8,-32)], hjust = justify_h[c(-23,-10,-8,-32)], 
	label = wrap_dACC[c(-23,-10,-8,-32)], lineheight = 0.85,color=text_color, family=font, vjust=0.5) +
annotate("text", y = all_dACC_sig2$Description[c(8,10)], x=c(-2.72,-2.12) , hjust = 1, #position_nudge(x=-0.1),
	label = wrap_dACC[c(8,10)], lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO dACC", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-10.5,2.5),labels = c(seq(10.5, 0, by = -0.5),seq(0.5, 2.5, by = 0.5)), breaks=seq(-10.5, 2.5, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.68),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(t=0.25,r=0,b=9,l=0),hjust = 0.68),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5), hjust=0.70),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold",margin = margin(r=-50)),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=dACC_dumbell,file="png/GO_dACC_AcrossDisorders_dumbbell_latest.png",width=10,height=8.5,dpi=300)


#############################################################################################
#DLPFC
#############################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
library(ggplot2)
library(RColorBrewer)
library(scales)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSD.rda")
go_DLPFC_PTSD <- go_DLPFC
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_onlyPTSD.rda")
go_DLPFC_onlyPTSD <- go_DLPFC
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_MDD.rda")
go_DLPFC_MDD <- go_DLPFC
load("rdas/DLPFC/geneSet_threeGroups_qSVA_DLPFC_PTSDvsMDD.rda")
go_DLPFC_PTSDvsMDD <- go_DLPFC

#change to data frame
go_DLPFC_PTSD_df <- as.data.frame(go_DLPFC_PTSD)
go_DLPFC_onlyPTSD_df <- as.data.frame(go_DLPFC_onlyPTSD)
go_DLPFC_MDD_df <- as.data.frame(go_DLPFC_MDD)
go_DLPFC_PTSDvsMDD_df <- as.data.frame(go_DLPFC_PTSDvsMDD)

#keep only significant (adjusted p < 0.1)
go_DLPFC_PTSD_sig <- go_DLPFC_PTSD_df[go_DLPFC_PTSD_df$p.adjust < 0.1,] 
go_DLPFC_onlyPTSD_sig <- go_DLPFC_onlyPTSD_df[go_DLPFC_onlyPTSD_df$p.adjust < 0.1,]
go_DLPFC_MDD_sig <- go_DLPFC_MDD_df[go_DLPFC_MDD_df$p.adjust < 0.1,]
go_DLPFC_PTSDvsMDD_sig <- go_DLPFC_PTSDvsMDD_df[go_DLPFC_PTSDvsMDD_df$p.adjust < 0.1,] #none sig

#uniquely identify the three analyses so can dstinguish after rbinding them
go_DLPFC_PTSD_sig$Analysis = "PTSD"
go_DLPFC_onlyPTSD_sig$Analysis = "onlyPTSD"
go_DLPFC_MDD_sig$Analysis = "MDD" 
#go_DLPFC_PTSDvsMDD_sig$Analysis = "PTSDvsMDD"

#add go_DLPFC_PTSDvsMDD_sig if any sig
all_DLPFC_sig <- rbind(go_DLPFC_onlyPTSD_sig, go_DLPFC_MDD_sig, go_DLPFC_PTSD_sig)

#still need to incorporate up vs down
#get rid of p005all for this graph
all_DLPFC_sig <- all_DLPFC_sig[all_DLPFC_sig$Cluster == "p005up" | all_DLPFC_sig$Cluster == "p005down",]
#make new column with -log10(adjP) calculated 
all_DLPFC_sig$P.Adj.Calc <- (-log10(all_DLPFC_sig$p.adjust))
#make another column with the calculated value as negative if p005 down, positive if p005up
all_DLPFC_sig$final.P.Adj <- ifelse(all_DLPFC_sig$Cluster == "p005up", all_DLPFC_sig$P.Adj.Calc, -all_DLPFC_sig$P.Adj.Calc)


all_DLPFC_sig[which(all_DLPFC_sig$Description %in% 
	intersect(all_DLPFC_sig$Description[all_DLPFC_sig$Analysis == "PTSD"],intersect(all_DLPFC_sig$Description[all_DLPFC_sig$Analysis=="onlyPTSD"],
	all_DLPFC_sig$Description[all_DLPFC_sig$Analysis=="MDD"]))),]

intersect(all_DLPFC_sig$Description[all_DLPFC_sig$Analysis == "PTSD"],all_DLPFC_sig$Description[all_DLPFC_sig$Analysis=="onlyPTSD"])

intersect(all_DLPFC_sig$Description[all_DLPFC_sig$Analysis == "MDD"],all_DLPFC_sig$Description[all_DLPFC_sig$Analysis=="PTSD"])                                         

all_DLPFC_sig$GeneRatioDec <- sapply(all_DLPFC_sig$GeneRatio, function(x) eval(parse(text = x)))
all_DLPFC_sig <- all_DLPFC_sig[order(all_DLPFC_sig$GeneRatioDec, decreasing=TRUE),]

all_DLPFC_sig[all_DLPFC_sig$Analysis == "onlyPTSD",][,"Description"]
#all_DLPFC_sig[all_DLPFC_sig$Analysis == "PTSDvsMDD",][1:20,"Description"]
all_DLPFC_sig[all_DLPFC_sig$Analysis == "PTSD",][,"Description"]
all_DLPFC_sig[all_DLPFC_sig$Analysis == "MDD",][,"Description"]

DLPFC_keep <- c("response to ketone","steroid hormone mediated signaling pathway","estrogen receptor binding","cellular response to lipid",
	"peroxidase activity","nuclear-transcribed mRNA catabolic process, nonsense-mediated decay","regulation of plasma lipoprotein particle levels",
	"ribosomal subunit","translational initiation","SRP-dependent cotranslational protein targeting to membrane","protein localization to endoplasmic reticulum",
	"serine-type endopeptidase activity","protein-lipid complex assembly","protein transport along microtubule","secretory granule localization",
	"establishment of protein localization to organelle","retrograde axonal transport")

all_DLPFC_sig2 <- all_DLPFC_sig[which(all_DLPFC_sig$Description %in% DLPFC_keep),]

#dumbbell

library("ggalt")
library("tidyr")
library(dplyr)
library(stringr)

x_colors<- c(rep("#C24A4A",4),"#222222",rep("#5A9C64",4))
font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("onlyPTSD" = "#3881B0","MDD" = "#CB8CAD", "PTSD" = "#000000")

all_DLPFC_sig2<- all_DLPFC_sig2[order(all_DLPFC_sig2$final.P.Adj,decreasing=FALSE),]

all_DLPFC_sig2$Description <- factor(all_DLPFC_sig2$Description, levels = unique(all_DLPFC_sig2$Description))

all_DLPFC_sig2$point_sizes<- ifelse(all_DLPFC_sig2$GeneRatioDec < 0.01, NA,ifelse(all_DLPFC_sig2$GeneRatioDec >= 0.01 & all_DLPFC_sig2$GeneRatioDec < 0.05, 4,
	 ifelse(all_DLPFC_sig2$GeneRatioDec >= 0.05 & all_DLPFC_sig2$GeneRatioDec < 0.1, 6, ifelse(all_DLPFC_sig2$GeneRatioDec >= 0.1 & all_DLPFC_sig2$GeneRatioDec < 0.15, 8,
	ifelse(all_DLPFC_sig2$GeneRatioDec >= 0.15 & all_DLPFC_sig2$GeneRatioDec < 0.2, 10, NA)))))

#save object so that i can make figure as png on my laptop (off the cluster)
save(all_DLPFC_sig2, file="rdas/GO_DLPFC_p.adjust0.1_plotinput.rda")

##The below generates a plot with the dumbbell end colored by Analysis and the segment is a solid grey line

pdf("pdf/GO_DLPFC_AcrossDisorders_dumbbell_latest.pdf",width=10)
ggplot(data=all_DLPFC_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_DLPFC_sig2,aes(y=all_DLPFC_sig2$Description, x=all_DLPFC_sig2$final.P.Adj,color=Analysis),size=4.25, position=position_nudge(y=0.05)) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_DLPFC_sig2$Description, x=ifelse(all_DLPFC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_DLPFC_sig2$final.P.Adj < 0, 0, 1), 
	label = stringr::str_wrap(all_DLPFC_sig2$Description,70),vjust=0.5, lineheight = 0.85,color=text_color, family=font) +
ggtitle("GO DLPFC", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2,2),labels = c(seq(2, 0, by = -0.5),seq(0.5, 2, by = 0.5)), breaks=seq(-2, 2, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(family=font,size=18,color=text_color,face="bold"),
	legend.title.align = 0.5,
	legend.text = element_text(family=font,size=16,color=text_color),
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
dev.off()

#png
#i have to generate this on my laptop, not the cluster, since ggsave has font issues with the cluster

x_colors<- c(rep("#C24A4A",4),"#222222",rep("#5A9C64",4))
font <- "Helvetica"
text_color <- "#222222"
analysis_colors <- c("onlyPTSD" = "#3881B0","MDD" = "#CB8CAD", "PTSD" = "#000000")
point_sizes<-all_DLPFC_sig2$point_sizes 	
DLPFC_wrap<-stringr::str_wrap(all_DLPFC_sig2$Description,70)
DLPFC_wrap[15] <- "nuclear-transcribed mRNA catabolic process,\nnonsense-mediated decay"

DLPFC_dumbell <- ggplot(data=all_DLPFC_sig2, aes(x=0, xend=final.P.Adj,y = Description)) +
geom_dumbbell(colour="#dddddd", colour_x = NULL, colour_xend=NULL, size_x = 0, size_xend = 0, size=2.5) + 
geom_point(data=all_DLPFC_sig2,aes(y=all_DLPFC_sig2$Description, x=all_DLPFC_sig2$final.P.Adj,color=Analysis),size=point_sizes) +
scale_colour_manual(values=analysis_colors) +
annotate("text", y = all_DLPFC_sig2$Description, x=ifelse(all_DLPFC_sig2$final.P.Adj < 0,0.1,-0.1), hjust = ifelse(all_DLPFC_sig2$final.P.Adj < 0, 0, 1), 
	label = DLPFC_wrap, lineheight = 0.85,color=text_color, family=font, vjust=0.5) +
ggtitle("GO DLPFC", subtitle="Across Analyses") + 
xlab("-log10(Adjusted P Value)") +
ylab("Description") +
scale_x_continuous(limits=c(-2,2),labels = c(seq(2, 0, by = -0.5),seq(0.5, 2, by = 0.5)), breaks=seq(-2, 2, by = 0.5)) +
theme(
	panel.background = element_blank(),
	plot.title = element_text(family=font,size=24,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=20,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=18,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=18,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=16,color=x_colors, face="bold"),
	axis.ticks.x=element_line(colour = x_colors, size=1.15),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
ggsave(plot=DLPFC_dumbell,file="png/GO_DLPFC_AcrossDisorders_dumbbell_latest.png",width=10,height=8.5,dpi=300)


