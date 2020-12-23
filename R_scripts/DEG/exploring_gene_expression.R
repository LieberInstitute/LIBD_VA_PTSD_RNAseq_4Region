library(jaffelab)
library(SummarizedExperiment)
library(ggplot2)
library(RColorBrewer)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

tog <- as.data.frame(assays(rse_gene)$counts[rownames(assays(rse_gene)$counts) == "ENSG00000241563.3"])
tog$Region <- rse_gene$Region
tog$Group <- rse_gene$Group
names(tog) <- c("counts","Region","Group")


cort<-ggplot(tog,aes(x=Group,y=log2(counts + 1),color=Region)) +
geom_boxplot() +
theme(
	legend.background = element_blank(),
	legend.box.background = element_blank(), 
	legend.key=element_blank(),
	legend.title = element_text(size=18,face="bold"),
	legend.text = element_text(size=16,face="bold"),
	legend.title.align = 0.5,
	axis.title.x = element_text(size=18,face="bold"),
	axis.title.y = element_text(size=18,face="bold"),
	axis.text.x = element_text(size=16,face="bold"),
	axis.ticks.x=element_line(size=1),
	axis.ticks.length.x=unit(0.25,"cm"),
	axis.ticks.length.y=unit(0.25,"cm"),
	axis.text.y=element_text(size=16,face="bold")
)
ggsave(plot=cort, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/png/CORT_counts_boxplot.png")
