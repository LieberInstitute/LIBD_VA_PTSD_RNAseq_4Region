#Manhattan plots for DEGs
#BasoAmyg

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')
library(recount)
library(limma)
library(edgeR)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')


#load rse object and get into form that was used for DE analysis
load('rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
rse_gene$onlyPTSD = ifelse(rse_gene$Group == "PTSD" ,1, 0) 

## add MDS (get ethnicity via genotype)
load("rdas/PTSD_LIBD_VA_MDSonly_n326.rda")
rownames(mds) = ss(rownames(mds),"_")

#gene
colData(rse_gene) = cbind(colData(rse_gene) , mds[rse_gene$BrNum,])

#gene
gIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.2
rse_gene <- rse_gene[gIndex , ]

#filter for BasoAmyg 
#gene
keepIndex = which(rse_gene$Region == "BasoAmyg")
rse_gene <- rse_gene[, keepIndex]

#extract phenotype data
pd = colData(rse_gene)
#extract rowRanges
rr = data.frame(rowRanges(rse_gene))

#load DE results
load("rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda")

#break DE results into MDD, PTSD, onlyPTSD, PTSDvsMDD, ANOVA
MDD <- data.frame(geneStats_BasoAmygall[,grep("_MDD",colnames(geneStats_BasoAmygall))])
PTSD <- data.frame(geneStats_BasoAmygall[,grep("_PTSD\\b",colnames(geneStats_BasoAmygall))])
PTSDvsMDD <- data.frame(geneStats_BasoAmygall[,grep("_PTSDvsMDD",colnames(geneStats_BasoAmygall))])
onlyPTSD <- data.frame(geneStats_BasoAmygall[,grep("_onlyPTSD",colnames(geneStats_BasoAmygall))])
ANOVA <- data.frame(geneStats_BasoAmygall[,grep("_ANOVA",colnames(geneStats_BasoAmygall))])

#make required objects for ManhattanPlotFunctions.R

#get dataframe of chromosome sizes with column names "chr" and "size"
#using Human Genome Assembly GRCh38.p7 from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38.p7
load("rdas/chrom_sizes_GRCh38P7.rda")

#make df of differential expression analysis result with 5 columns
#    - 1st column : unique gene identifier (gene_id, external_gene_name etc.) --> ensemblID
#    - 2nd column : the probability measure of significance (P.Value, adj.P.Val etc) --> adj.P.Val
#    - 3rd/4th/5th column : to be named "chr", "start" and "end" respectively.    
#    - Note: column names of first two columns not important

length(unique(PTSD$ensemblID_PTSD))
length(PTSD$ensemblID_PTSD)
rownames(rr) <- rr$ensemblID
rownames(MDD) <- MDD$ensemblID

#test <- cbind(MDD, rr[MDD$ensemblID_MDD,])
#identical(test$ensemblID_MDD, test$ensemblID)
MDD_app <- cbind(MDD, rr[MDD$ensemblID_MDD,c(1:3)])
PTSD_app <- cbind(PTSD, rr[PTSD$ensemblID_PTSD,c(1:3)])
PTSDvsMDD_app <- cbind(PTSDvsMDD, rr[PTSDvsMDD$ensemblID_PTSDvsMDD,c(1:3)])
onlyPTSD_app <- cbind(onlyPTSD, rr[onlyPTSD$ensemblID_onlyPTSD,c(1:3)])
ANOVA_app <- cbind(ANOVA, rr[ANOVA$ensemblID_ANOVA,c(1:3)])

#get rid of extra stuff
MDD_short <- MDD_app[,c(3,15,17,18,19)]
PTSD_short <- PTSD_app[,c(3,15,17,18,19)]
PTSDvsMDD_short <- PTSDvsMDD_app[,c(3,15,17,18,19)]
onlyPTSD_short <- onlyPTSD_app[,c(3,15,17,18,19)]
ANOVA_short <- ANOVA_app[,c(3,15,17,18,19)]

#rename seqnames to chr
colnames(MDD_short)[3] <- "chr"
colnames(PTSD_short)[3] <- "chr"
colnames(PTSDvsMDD_short)[3] <- "chr"
colnames(onlyPTSD_short)[3] <- "chr"
colnames(ANOVA_short)[3] <- "chr"

# 3) bin size in base pairs to be used (default 1000000)
# 4) the probability cutoff for determining significance (default 0.05)
# 5) prefix for output files (default "manplot")

colnames(chrom_sizes)[1] <- "chr"
colnames(chrom_sizes)[2] <- "size"
MDD_short$chr=sub("chr","",MDD_short$chr)

chrsizetab <- chrom_sizes
detab <- MDD_short

plotManhattanPlotRNAseq<-function(chrsizetab,detab,binsize=1000000,pcutoff=0.5,prefix="manplot") {
# order detab by probability, most significant at the top
detab=detab[order(detab[,2],decreasing=F),]
head(detab)
# check for duplicated gene identifiers and remove the less significant one
table(duplicated(detab[,1]))
detab=detab[!duplicated(detab[,1]),]
# set gene midpoint, to be used as the representative position of the gene
detab$mid=(detab$start+detab$end)/2
# remove chrM from being considered
detab=detab[grep("chrM",detab$chr,invert=T),]
chrsizetab=chrsizetab[grep("chrM",chrsizetab$chr,invert=T),]
# more pre-processing of chrsizetab
# ensure size column is numeric 
chrsizetab$size=as.numeric(chrsizetab$size) 
# chrnum is the numeric chromosome number so that it can be sorted
chrsizetab$chrnum=sub("chr","",chrsizetab$chr)
chrsizetab$chrnum=sub("Y",length(row.names(chrsizetab)),chrsizetab$chrnum)
chrsizetab$chrnum=sub("X",length(row.names(chrsizetab))-1,chrsizetab$chrnum)
chrsizetab$chrnum=as.numeric(chrsizetab$chrnum)
chrsizetab=chrsizetab[order(chrsizetab$chrnum,decreasing=F),]
# calculate number of bins per chromosome
chrsizetab$numbins=ceiling(chrsizetab$size/binsize)
# cumnumbins is the cumulative number of bins
chrsizetab[1,"cumnumbins"]=chrsizetab[1,"numbins"]
# xpos is the middle bin of each chromosome to place the x axis labels
chrsizetab[1,"xpos"]=chrsizetab[1,"numbins"]/2
for (x in 2:length(row.names(chrsizetab))) {
chrsizetab[x,"cumnumbins"]=chrsizetab[x,"numbins"]+chrsizetab[x-1,"cumnumbins"]
chrsizetab[x,"xpos"]=(chrsizetab[x,"numbins"]/2)+(chrsizetab[x-1,"cumnumbins"])
}
chrsizetab
degtab=detab[detab[,2]<=pcutoff,] 
nondegtab=detab[detab[,2]>pcutoff,]
head(degtab)
head(nondegtab)
# assign each gene to a bin based on position
degtab$bin=paste(degtab$chr,ceiling(degtab$mid/binsize),sep="_")
head(degtab)
# create a data frame counting the number of DEG per bin
degdf=data.frame(table(degtab$bin))
names(degdf)=c("bin","deFreq")
head(degdf)
# store the DEG names as well
degdf$degenes=NA
for (row in 1:length(row.names(degtab))) {
#row=1
if (is.na(degdf[degdf$bin==degtab$bin[row],"degenes"])) {
degdf[degdf$bin==degtab$bin[row],"degenes"]=degtab[row,1]
}else {
degdf[degdf$bin==degtab$bin[row],"degenes"]=paste(degdf[degdf$bin==degtab$bin[row],"degenes"],degtab[row,1],sep="::")
    }
  }
# similar for nondegtab, assign each gene to a bin based on position, and create a data 
# frame counting the number of non-DEG per bin
nondegtab$bin=paste(nondegtab$chr,ceiling(nondegtab$mid/binsize),sep="_")
head(nondegtab)
nondegdf=data.frame(table(nondegtab$bin))
names(nondegdf)=c("bin","nondeFreq") 
head(nondegdf)
# merge degdf and nondegdf
tempdf=merge(degdf,nondegdf,by.x="bin",by.y="bin",all=T)  
tempdf[is.na(tempdf)]=0
tempdf[tempdf$degenes==0,"degenes"]=NA
head(tempdf)
for (j in 1:length(row.names(tempdf))) {
#j=1
tempdf[j,"deFreq"]
tempdf[j,"nondeFreq"]
tempdf[j,"deFreq_outsidebin"]=sum(tempdf[-j,"deFreq"])
tempdf[j,"nondeFreq_outsidebin"]=sum(tempdf[-j,"nondeFreq"])
m=matrix(c(tempdf[j,"deFreq"],tempdf[j,"nondeFreq"],tempdf[j,"deFreq_outsidebin"],tempdf[j,"nondeFreq_outsidebin"]),nrow=2,byrow=T)
f=fisher.test(m)
tempdf[j,"p.value"]=f$p.value
#cat(j," ",sep="")
  }
tempdf$p.adj=p.adjust(tempdf$p.value,method="BH")
head(tempdf)
#sum(tempdf$p.value<0.05)
#sum(tempdf$p.adj<0.05)
# create and populate a plotdf data frame with a row for every bin, ordered from chr1
# to chrY
chrsizetab
plotdf=data.frame()
for (i in 1:length(row.names(chrsizetab))) {
for (j in 1:chrsizetab[i,"numbins"]) {
plotdf[paste(chrsizetab[i,"chr"],j,sep="_"),"chr"]=sub("chr","",chrsizetab[i,"chr"])
plotdf[paste(chrsizetab[i,"chr"],j,sep="_"),"binnum"]=j
    }
  }
plotdf$chrname=sub("\\_\\d+","",row.names(plotdf))
# merge the various gene counts for bins that contain genes 
plotdf=merge(plotdf,tempdf,by.x=0,by.y="bin",all.x=T)
# comment out the following two lines if you do not want bins that are empty of genes to
# have a point on the manhattan plot.
#plotdf[is.na(plotdf$p.value),"p.value"]=1
#plotdf[is.na(plotdf$p.adj),"p.adj"]=1
plotdf$neglog10pval= -log(plotdf$p.value,10)
plotdf$neglog10padj= -log(plotdf$p.adj,10)
plotdf[plotdf$chr=="X","chr"]=length(row.names(chrsizetab))-1
plotdf[plotdf$chr=="Y","chr"]=length(row.names(chrsizetab))
plotdf$chr=as.numeric(plotdf$chr)
plotdf=plotdf[order(plotdf$chr,plotdf$binnum,decreasing=F),]
row.names(plotdf)=1:length(row.names(plotdf))
plotdf$color=rainbow(length(row.names(chrsizetab)))[plotdf$chr]
plotdf$binstartpos=(plotdf$binnum-1)*binsize
plotdf$binendpos=(plotdf$binnum)*binsize
plotdf$binmidpos=(plotdf$binstartpos+plotdf$binendpos)/2
head(plotdf)
summary(plotdf$neglog10pval)
summary(plotdf$neglog10padj)
write.table(plotdf,paste(prefix,"_",(binsize/1000000),"Mb-bin_data.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
pdf(paste(prefix,"_",(binsize/1000000),"Mb-bin.pdf",sep=""),width=11,height=8.5, onefile=TRUE)
par(mfrow = c(2, 1),oma=c(1,1,1,1),mar=c(4,4,4,4)+0.1)
ymax=ceiling(max(c(2,plotdf$neglog10pval),na.rm=T))
plot(1:length(row.names(plotdf)),plotdf$neglog10pval,ylim=c(0,ymax),col=plotdf$color,bty="n",ylab="-log10(p.value)",xlab="Chromosome",xaxt="n",main=paste(prefix," Manhattan Plot\n(binsize=",(binsize/1000000)," Mb)",sep=""))
abline(h= (-log(0.05,10)),lty=2)
axis(side=1,at=chrsizetab$xpos,labels=sub("chr","",chrsizetab$chr))
ymax=ceiling(max(c(2,plotdf$neglog10padj),na.rm=T))
plot(1:length(row.names(plotdf)),plotdf$neglog10padj,ylim=c(0,ymax),col=plotdf$color,bty="n",ylab="-log10(p.adj)",xlab="Chromosome",xaxt="n",main=paste(prefix," Manhattan Plot\n(binsize=",(binsize/1000000)," Mb)",sep=""))
abline(h= (-log(0.05,10)),lty=2)
axis(side=1,at=chrsizetab$xpos,labels=sub("chr","",chrsizetab$chr))
par(mfrow = c(1, 1),oma=c(1,1,1,1),mar=c(4,4,4,4)+0.1)
graphics.off()
}




