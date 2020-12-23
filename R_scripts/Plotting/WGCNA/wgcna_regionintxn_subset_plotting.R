
#Boxplots of eigengene vs dx for modules that significantly associated with dx

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(devtools)
library(cluster)
library(clusterProfiler)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

########################################################################
####Get diagnosis for each RNum
########################################################################

#load rse objects
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

rse_gene$CatRegion <- ifelse(rse_gene$Region == "BasoAmyg" | rse_gene$Region == "MedialAmyg", "Amyg","Cortex")

pd<-colData(rse_gene)

###Amyg

#MDD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_MDD.rda",verbose=TRUE)
Amyg_MDD_MEs <- net$MEs
colnames(Amyg_MDD_MEs) <- paste0("MDD_",colnames(Amyg_MDD_MEs))

#PTSD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_PTSD.rda",verbose=TRUE)
Amyg_PTSD_MEs <- net$MEs
colnames(Amyg_PTSD_MEs) <- paste0("PTSD_",colnames(Amyg_PTSD_MEs))

#onlyPTSD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Amyg_onlyPTSD.rda",verbose=TRUE)
Amyg_onlyPTSD_MEs <- net$MEs
colnames(Amyg_onlyPTSD_MEs) <- paste0("onlyPTSD_",colnames(Amyg_onlyPTSD_MEs))

Amyg_MEs <- cbind(Amyg_MDD_MEs, Amyg_PTSD_MEs, Amyg_onlyPTSD_MEs)
Amyg_MEs$Dx <- pd$Group[match(rownames(Amyg_MEs),rownames(pd))]


###Cortex

#MDD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_MDD.rda",verbose=TRUE)
Cortex_MDD_MEs <- net$MEs
colnames(Cortex_MDD_MEs) <- paste0("MDD_",colnames(Cortex_MDD_MEs))

#PTSD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_PTSD.rda",verbose=TRUE)
Cortex_PTSD_MEs <- net$MEs
colnames(Cortex_PTSD_MEs) <- paste0("PTSD_",colnames(Cortex_PTSD_MEs))

#onlyPTSD modules
load("rdas/WGCNA/constructed_network_signed_bicor_Cortex_onlyPTSD.rda",verbose=TRUE)
Cortex_onlyPTSD_MEs <- net$MEs
colnames(Cortex_onlyPTSD_MEs) <- paste0("onlyPTSD_",colnames(Cortex_onlyPTSD_MEs))

Cortex_MEs <- cbind(Cortex_MDD_MEs, Cortex_PTSD_MEs, Cortex_onlyPTSD_MEs)
Cortex_MEs$Dx <- pd$Group[match(rownames(Cortex_MEs),rownames(pd))]

save(Amyg_MEs,Cortex_MEs,file="rdas/WGCNA/WGCNA_regionintxn_subset_MEs_for_plotting.rda")


##############################################
###Check out sig module eigengenes
##############################################
#Amyg
load("rdas/WGCNA/MEvsDx_Amyg_MDD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME20
rownames(coefAdj) <- paste0("MDD_",rownames(coefAdj))
sig_amyg_MDD <- coefAdj[coefAdj[,5] < 0.05,5]

load("rdas/WGCNA/MEvsDx_Amyg_PTSD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME2
rownames(coefAdj) <- paste0("PTSD_",rownames(coefAdj))
sig_amyg_PTSD <- coefAdj[coefAdj[,5] < 0.05,5]

load("rdas/WGCNA/MEvsDx_Amyg_onlyPTSD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME2
rownames(coefAdj) <- paste0("onlyPTSD_",rownames(coefAdj))
sig_amyg_onlyPTSD <- coefAdj[coefAdj[,5] < 0.05,5]
names(sig_amyg_onlyPTSD) <- "onlyPTSD_ME2"

sig_amyg <- c(sig_amyg_MDD, sig_amyg_PTSD, sig_amyg_onlyPTSD)


#Cortex
load("rdas/WGCNA/MEvsDx_Cortex_MDD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME8
rownames(coefAdj) <- paste0("MDD_",rownames(coefAdj))
sig_cortex_MDD <- coefAdj[coefAdj[,5] < 0.05,5]

load("rdas/WGCNA/MEvsDx_Cortex_PTSD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME10
rownames(coefAdj) <- paste0("PTSD_",rownames(coefAdj))
sig_cortex_PTSD <- coefAdj[coefAdj[,5] < 0.05,5]

load("rdas/WGCNA/MEvsDx_Cortex_onlyPTSD.rda",verbose=TRUE)
coefAdj[which(coefAdj[,5] %in% min(coefAdj[,5])),,drop=FALSE]
#ME3
rownames(coefAdj) <- paste0("onlyPTSD_",rownames(coefAdj))
sig_cortex_onlyPTSD <- coefAdj[coefAdj[,5] < 0.05,5]

sig_cortex <- c(sig_cortex_MDD, sig_cortex_PTSD, sig_cortex_onlyPTSD)

save(sig_amyg,sig_cortex,file="rdas/WGCNA/WGCNA_regionintxn_subset_sig_ME_names_for_plotting.rda")

######################################
#Get GO terms
######################################

wgcna_MEs_GO <- vector(mode="list",length =2)
names(wgcna_MEs_GO)<-c("Amyg","Cortex")

wgcna_MEs_GO_sig <- vector(mode="list",length =2)
names(wgcna_MEs_GO_sig)<-c("Amyg","Cortex")

##Amyg

#MDD
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_MDD.rda",verbose=TRUE)
go_modules_amyg_MDD <- as.data.frame(go_modules_amyg)
go_modules_amyg_MDD$Cluster <- paste0("MDD_ME",go_modules_amyg_MDD$Cluster)

#PTSD
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_PTSD.rda",verbose=TRUE)
go_modules_amyg_PTSD <- as.data.frame(go_modules_amyg)
go_modules_amyg_PTSD$Cluster <- paste0("PTSD_ME",go_modules_amyg_PTSD$Cluster)

#onlyPTSD
load("rdas/WGCNA/wgcna_Amyg_GO_clusterProfiler_onlyPTSD.rda",verbose=TRUE)
go_modules_amyg_onlyPTSD <- as.data.frame(go_modules_amyg)
go_modules_amyg_onlyPTSD$Cluster <- paste0("onlyPTSD_ME",go_modules_amyg_onlyPTSD$Cluster)

go_modules_amyg_all<-rbind(go_modules_amyg_MDD, go_modules_amyg_PTSD, go_modules_amyg_onlyPTSD)
wgcna_MEs_GO$Amyg <- split(go_modules_amyg_all,go_modules_amyg_all$Cluster)

wgcna_MEs_GO_sig$Amyg <- wgcna_MEs_GO$Amyg[which(names(wgcna_MEs_GO$Amyg) %in% names(sig_amyg))]

##Cortex

#MDD
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_MDD.rda",verbose=TRUE)
go_modules_cortex_MDD <- as.data.frame(go_modules_cortex)
go_modules_cortex_MDD$Cluster <- paste0("MDD_ME",go_modules_cortex_MDD$Cluster)

#PTSD
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_PTSD.rda",verbose=TRUE)
go_modules_cortex_PTSD <- as.data.frame(go_modules_cortex)
go_modules_cortex_PTSD$Cluster <- paste0("PTSD_ME",go_modules_cortex_PTSD$Cluster)

#onlyPTSD
load("rdas/WGCNA/wgcna_Cortex_GO_clusterProfiler_onlyPTSD.rda",verbose=TRUE)
go_modules_cortex_onlyPTSD <- as.data.frame(go_modules_cortex)
go_modules_cortex_onlyPTSD$Cluster <- paste0("onlyPTSD_ME",go_modules_cortex_onlyPTSD$Cluster)

go_modules_cortex_all<-rbind(go_modules_cortex_MDD, go_modules_cortex_PTSD, go_modules_cortex_onlyPTSD)
wgcna_MEs_GO$Cortex <- split(go_modules_cortex_all,go_modules_cortex_all$Cluster)

wgcna_MEs_GO_sig$Cortex <- wgcna_MEs_GO$Cortex[which(names(wgcna_MEs_GO$Cortex) %in% names(sig_cortex))]

save(wgcna_MEs_GO,wgcna_MEs_GO_sig, file="rdas/WGCNA/WGCNA_regionintxn_subset_GO_for_plotting.rda")


####################################
##Let's make some plots
####################################

#Get GO terms ready for auto plotting (mapply)

amyg_names <- names(wgcna_MEs_GO_sig$Amyg)

amyg_GO_topdesc <- lapply(amyg_names, function(x) {
	wgcna_MEs_GO_sig$Amyg[[x]][which(wgcna_MEs_GO_sig$Amyg[[x]]$pvalue == min(wgcna_MEs_GO_sig$Amyg[[x]]$pvalue)),"Description"][1]
})

names(amyg_GO_topdesc) <-amyg_names
amyg_GO_topdesc<-unlist(amyg_GO_topdesc)


cortex_names <- names(wgcna_MEs_GO_sig$Cortex)

cortex_GO_topdesc <- lapply(cortex_names, function(x) {
	wgcna_MEs_GO_sig$Cortex[[x]][which(wgcna_MEs_GO_sig$Cortex[[x]]$pvalue == min(wgcna_MEs_GO_sig$Cortex[[x]]$pvalue)),"Description"][1]
})

names(cortex_GO_topdesc) <-cortex_names
cortex_GO_topdesc<-unlist(cortex_GO_topdesc)

#####Check out all sig

#Amygdala
Amyg_MEs <- Amyg_MEs[,c(which(colnames(Amyg_MEs) %in% names(sig_amyg)),73)]

Amyg_MEs_Dx<-data.frame()
Amyg_MEs_Dx[1:644,1:9]<-Amyg_MEs$Dx
names(Amyg_MEs_Dx) <- names(Amyg_MEs)[1:9]
Amyg_MEs_Dx_list <- as.list(Amyg_MEs_Dx)

Amyg_MEs_list <- Amyg_MEs[,-10]
Amyg_MEs_list <- as.list(Amyg_MEs_list)

amyg_GO<-amyg_GO_topdesc[match(names(Amyg_MEs_list), names(amyg_GO_topdesc))]

pdf('Figures/PDFs_other/WGCNA_MEs_Dx_boxplots_Amyg.pdf', useDingbats = FALSE)
mapply(function(x, y, module, region, go, pval) {
    dx <- y
	ylim <- range(Amyg_MEs_list)
    dx <- factor(ifelse(dx == 'PTSD', 'PTSD', ifelse(dx == 'MDD', 'MDD', 'Control')))
    boxplot(x ~ dx, main = paste(region, '-', module, "\nP Value =",signif(pval, 3), "\nTop GO:", go),
        xlab = 'Diagnosis', outline = FALSE, ylab = 'Module Eigengene', ylim=ylim)
    points(x ~ jitter(as.numeric(dx), amount = 0.2), cex = 1.5, pch = 21, bg = '#FFCF00')
}, Amyg_MEs_list, Amyg_MEs_Dx_list, names(Amyg_MEs_list), "Amygdala",amyg_GO,sig_amyg)
dev.off()


#Cortex
Cortex_MEs <- Cortex_MEs[,c(which(colnames(Cortex_MEs) %in% names(sig_cortex)),92)]

Cortex_MEs_Dx<-data.frame()
Cortex_MEs_Dx[1:641,1:29]<-Cortex_MEs$Dx
names(Cortex_MEs_Dx) <- names(Cortex_MEs)[1:29]
Cortex_MEs_Dx_list <- as.list(Cortex_MEs_Dx)

Cortex_MEs_list <- Cortex_MEs[,-30]
Cortex_MEs_list <- as.list(Cortex_MEs_list)

cortex_GO<-cortex_GO_topdesc[match(names(Cortex_MEs_list), names(cortex_GO_topdesc))]

pdf('Figures/PDFs_other/WGCNA_MEs_Dx_boxplots_Cortex.pdf', useDingbats = FALSE)
mapply(function(x, y, module, region, go, pval) {
    dx <- y
	ylim <- range(Cortex_MEs_list)
    dx <- factor(ifelse(dx == 'PTSD', 'PTSD', ifelse(dx == 'MDD', 'MDD', 'Control')))
    boxplot(x ~ dx, main = paste(region, '-', module, "\nP Value =",signif(pval, 3), "\nTop GO:", go),
        xlab = 'Diagnosis', outline = FALSE, ylab = 'Module Eigengene',ylim=ylim)
    points(x ~ jitter(as.numeric(dx), amount = 0.2), cex = 1.5, pch = 21, bg = '#C43EB4')
}, Cortex_MEs_list, Cortex_MEs_Dx_list, names(Cortex_MEs_list), "Cortex",cortex_GO,sig_cortex)
dev.off()


###Nice plots for modules of interest

load("rdas/WGCNA/WGCNA_regionintxn_subset_MEs_for_plotting.rda",verbose=TRUE)
load("rdas/WGCNA/WGCNA_regionintxn_subset_GO_for_plotting.rda",verbose=TRUE)

amyg_names <- names(wgcna_MEs_GO_sig$Amyg)
amyg_GO_topdesc <- lapply(amyg_names, function(x) {
	top<-wgcna_MEs_GO_sig$Amyg[[x]][order(wgcna_MEs_GO_sig$Amyg[[x]]$pvalue,decreasing=FALSE),"Description"]
	top[1:5]
})
names(amyg_GO_topdesc) <-amyg_names

cortex_names <- names(wgcna_MEs_GO_sig$Cortex)
cortex_GO_topdesc <- lapply(cortex_names, function(x) {
	top<-wgcna_MEs_GO_sig$Cortex[[x]][order(wgcna_MEs_GO_sig$Cortex[[x]]$pvalue,decreasing=FALSE),"Description"]
	top[1:5]
})
names(cortex_GO_topdesc) <-cortex_names



##Amyg: #FFCF00 -> #FFE8A6
#MDD modules: 20,22
#PTSD modules: 20,18,2

##Cortex: #C43EB4 -> #DAB0D3
#MDD modules: 5, 10, 8, 22 (very interesting...)
#PTSD modules: 1,10,7,8


library(ggplot2)
library(RColorBrewer)

#onlyPTSD: #3881B0
#MDD: #CB8CAD
#PTSD: #000000


text_color<-"#222222"
font<-"Helvetica"

##Amyg

#MDD, ME20
amyg_MDD20<-ggplot(data=Amyg_MEs,aes(x=Dx,y=MDD_ME20)) +
geom_boxplot(colour="#222222",fill ="#FFE8A6", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Heparan Sulfate Proteoglycan Metabolic Process, ECM Constituent Secretion, and Neuron Remodeling") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=amyg_MDD20, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Amyg_MDD_ME20.png",width=18,height=14,dpi=320)

#############MDD, ME22
amyg_MDD22<-ggplot(data=Amyg_MEs,aes(x=Dx,y=MDD_ME22)) +
geom_boxplot(colour="#222222",fill ="#FFE8A6", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
scale_y_continuous(limits=c(-0.25,0.2),labels = c(seq(-0.25, 0, by = 0.05),seq(0.05, 0.2, by = 0.05)), breaks=seq(-0.25, 0.2, by = 0.05)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Regulation of transmembrane transport and integral component of synaptic vesicle membrane") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=amyg_MDD22, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Amyg_MDD_ME22.png",width=18,height=14,dpi=320)


#PTSD, ME20
amyg_PTSD20<-ggplot(data=Amyg_MEs,aes(x=Dx,y=PTSD_ME20)) +
geom_boxplot(colour="#222222",fill ="#FFE8A6", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Regulation of system process, synapse structure or activity, and synapse assembly") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=amyg_PTSD20, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Amyg_PTSD_ME20.png",width=18,height=14,dpi=320)

##############PTSD, ME18
amyg_PTSD18<-ggplot(data=Amyg_MEs,aes(x=Dx,y=PTSD_ME18)) +
geom_boxplot(colour="#222222",fill ="#FFE8A6", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
scale_y_continuous(limits=c(-0.25,0.2),labels = c(seq(-0.25, 0, by = 0.05),seq(0.05, 0.2, by = 0.05)), breaks=seq(-0.25, 0.2, by = 0.05)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Extracellular matrix constituent secretion, neuron remodeling, and heparan sulfate sulfotransferase activity") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=amyg_PTSD18, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Amyg_PTSD_ME18.png",width=18,height=14,dpi=320)

#PTSD, ME2
amyg_PTSD2<-ggplot(data=Amyg_MEs,aes(x=Dx,y=PTSD_ME2)) +
geom_boxplot(colour="#222222",fill ="#FFE8A6", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Regulation of cell activation, lymphocyte activation, and activation of immune response") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=amyg_PTSD2, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Amyg_PTSD_ME2.png",width=18,height=14,dpi=320)


##Cortex: #C43EB4 -> #CD74C0
#MDD modules: 5, 10, 8, 22 (very interesting...)
#PTSD modules: 1,10,7,8

##Cortex 

#MDD, ME5
cortex_MDD5<-ggplot(data=Cortex_MEs,aes(x=Dx,y=MDD_ME5)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Rac protein signal transduction, behavior, and regulation of transmembrane transport") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_MDD5, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_MDD_ME5.png",width=18,height=14,dpi=320)

#MDD, ME10
cortex_MDD10<-ggplot(data=Cortex_MEs,aes(x=Dx,y=MDD_ME10)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Presynapse, glutamatergic synapse, and synaptic membrane") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_MDD10, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_MDD_ME10.png",width=18,height=14,dpi=320)

#MDD, ME8
cortex_MDD8<-ggplot(data=Cortex_MEs,aes(x=Dx,y=MDD_ME8)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Regulation of cell activation, activation of immune response, and lymphocyte activation") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_MDD8, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_MDD_ME8.png",width=18,height=14,dpi=320)


#########################MDD, ME22
cortex_MDD22<-ggplot(data=Cortex_MEs,aes(x=Dx,y=MDD_ME22)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#CB8CAD",size=6,position=position_jitter(0.2)) +
scale_y_continuous(limits=c(-0.25,0.2),labels = c(seq(-0.25, 0, by = 0.05),seq(0.05, 0.2, by = 0.05)), breaks=seq(-0.25, 0.2, by = 0.05)) +
ggtitle("WGCNA Module Associating with MDD", subtitle="Cerebral cortex GABAergic interneuron differentiation, neuropeptide signaling pathway, and forebrain neuron differentiation") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=22,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_MDD22, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_MDD_ME22.png",width=18,height=14,dpi=320)

###########################PTSD, ME1
cortex_PTSD1<-ggplot(data=Cortex_MEs,aes(x=Dx,y=PTSD_ME1)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
scale_y_continuous(limits=c(-0.25,0.2),labels = c(seq(-0.25, 0, by = 0.05),seq(0.05, 0.2, by = 0.05)), breaks=seq(-0.25, 0.2, by = 0.05)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Myelination, ensheathment of neurons, and axon ensheathment") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_PTSD1, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_PTSD_ME1.png",width=18,height=14,dpi=320)


#PTSD, ME10
cortex_PTSD10<-ggplot(data=Cortex_MEs,aes(x=Dx,y=PTSD_ME10)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Catalytic activity acting on RNA, nucleic acid phosphodiester bond hydrolysis, and tRNA metabolic process") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_PTSD10, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_PTSD_ME10.png",width=18,height=14,dpi=320)

#PTSD, ME7
cortex_PTSD7<-ggplot(data=Cortex_MEs,aes(x=Dx,y=PTSD_ME7)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="tRNA metabolic process, ncRNA processing, and catalytic activity acting on RNA") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_PTSD7, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_PTSD_ME7.png",width=18,height=14,dpi=320)

#PTSD, ME8
cortex_PTSD8<-ggplot(data=Cortex_MEs,aes(x=Dx,y=PTSD_ME8)) +
geom_boxplot(colour="#222222",fill ="#CD74C0", outlier.shape=NA,size=1.25) +
geom_jitter(shape=21,colour="#222222",fill="#000000",size=6,position=position_jitter(0.2)) +
ggtitle("WGCNA Module Associating with PTSD", subtitle="Regulation of cell activation, activation of immune response, and regulation of leukocyte activation") +
xlab("Diagnosis") +
ylab("Module Eigengene") + 
theme(
	panel.background = element_blank(),
	#axis.line=element_blank(),
	panel.grid.major.y = element_line(color="#D3D3D3",size=1),
	panel.grid.major.x = element_line(color="#D3D3D3",size=1),
    #panel.grid.minor.x = element_blank(),
	#panel.grid.minor.y = element_line(color="#222222",size=0.5),
    panel.border = element_rect(color="#222222",size=1.5,fill=NA),
	plot.title = element_text(family=font,size=28,face="bold",color=text_color,hjust = 0.5),
	plot.subtitle = element_text(family=font,size=24,margin=ggplot2::margin(0.25,0,9,0),hjust = 0.5),
	legend.position = "none",
	axis.title.x = element_text(family=font, size=22,color=text_color,face="bold", margin = margin(t = 10, b = 5)),
	axis.title.y = element_text(family=font, size=22,color=text_color,face="bold"),
	axis.text.x = element_text(family=font, size=20,color=text_color, face="bold"),
	axis.ticks.y=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.y=unit(0.5,"cm"),
	axis.ticks.x=element_line(colour="#D3D3D3",size=1),
	axis.ticks.length.x=unit(0.5,"cm"),
	axis.text.y = element_text(family=font, size=20,color=text_color, face="bold"))
ggsave(plot=cortex_PTSD8, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Figures/PNGs_other/WGCNA_MEs_Dx_boxplot_Cortex_PTSD_ME8.png",width=18,height=14,dpi=320)





