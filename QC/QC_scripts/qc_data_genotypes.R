###
library(rtracklayer)
library(SummarizedExperiment)
library(VariantAnnotation)
library(jaffelab)
library(readxl)

## load gene counts
load("preprocessed_data/rse_gene_PTSD_VA_LIBD_Combined_n1304.Rdata")
rse_gene$Region[rse_gene$Region == "Dacc"] = "dACC"
rse_gene$Region[rse_gene$Region == "BasolateralAmyg"] = "BasoAmyg"

## read in VCF
vcf = readVcf("preprocessed_data/mergedVariants_alreadyCombined.vcf.gz")
vcf = vcf[,rse_gene$bamFile] # match up to correct BAM file
colnames(vcf) = rse_gene$RNum

# add rs num
snpMap = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")
oo = findOverlaps(vcf, snpMap, type="equal")
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)]
rowData(vcf)$snpRsNum[is.na(rowData(vcf)$snpRsNum)] = rownames(vcf)[is.na(rowData(vcf)$snpRsNum)]
rownames(vcf) = rowData(vcf)$snpRsNum

#########################
## line up to pcr data
geno = as.data.frame(read_excel("VA_20brains_20180908.xlsx"))
geno = DataFrame(geno)
geno$Genotypes = strsplit(geno$SNPmark, "\r\n")
geno$Genotypes = lapply(geno$Genotypes, function(x) {
	y = ss(x, ": ", 2)
	names(y) = ss(x, ": ", 1)
	y
})

## take the subset in common
tt = table(unlist(lapply(geno$Genotypes, names)))
keepSnps = names( tt[tt==80])

genoMat = do.call("cbind", lapply(geno$Genotypes,
	function(x) x[keepSnps])) 
colnames(genoMat) = geno$RNum
genoMat[genoMat %in% c("NoCall", "Confilict")] = NA

## get same SNPs from RNA
vcfMatch = vcf[rownames(genoMat), colnames(genoMat)]

## matched SNPs
snpMatch = geno(vcfMatch)$GT

## flag
flagInd = colnames(vcfMatch) %in% c("R16766", "R16763", "R16309", "R16418", "R16436", 
	"R16437", "R16487", "R16437", "R16487","R16100", "R16776", 
	"R17026", "R16253", "R16275", "R17284", "R16558")

table(snpMatch[1,], genoMat[1,])
table(snpMatch[1,!flagInd], genoMat[1,!flagInd])
table(snpMatch[1,flagInd], genoMat[1,flagInd])

genoCode = data.frame(paste0(ref(vcfMatch),"/", ref(vcfMatch)),
	paste0(unlist(alt(vcfMatch)),"/", ref(vcfMatch)),
	paste0(unlist(alt(vcfMatch)),"/", unlist(alt(vcfMatch))),
	stringsAsFactors=FALSE)
snpMatchA = snpMatch
for(i in 1:nrow(snpMatchA)) {
	xx= snpMatchA[i,] 
	xx[xx=="."] = genoCode[i,1]
	xx[xx=="0/1"] = genoCode[i,2]
	xx[xx=="1/1"] = genoCode[i,3]
	snpMatchA[i,] = xx
}

## check
numCheck = vector("numeric", ncol(snpMatchA))
for(i in seq(along=numCheck)) {
	numCheck[i] = table(snpMatchA[,i] ==  genoMat[,i], useNA="ifany")[2]
}
boxplot(numCheck ~  flagInd)
numCheck[which(lengths(apply(genoMat, 1, table)) == 3)]

table(snpMatchA[1,] == genoMat[1,])

table(snpMatchA[1,], genoMat[1,])
table(snpMatchA[5,], genoMat[5,])
table(snpMatchA[1,!flagInd], genoMat[1,!flagInd])
table(snpMatchA[1,flagInd], genoMat[1,flagInd])

## just observed
genoMatFac = as.data.frame(genoMat)
genoMatNum = genoMatFac
for(i in 1:ncol(genoMatNum)) genoMatNum[,i] = as.numeric(genoMatNum[,i])
colnames(genoMatNum) = paste0(colnames(genoMatNum),":", geno$BrNum[match(colnames(genoMatNum), geno$RNum)])

ccObs = cor(genoMatNum, use="comp")
hcObs = hclust(as.dist(1-ccObs))
myplclust(hcObs)

#########################
	######################
# subset to high-depth
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          info(vcf)$VDB >0.1,]
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
		 
########################################
# check snp correlation of all samples
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snpCor = cor(snps, use="pairwise.complete.obs")

### PLOT ####
library(pheatmap)
library(RColorBrewer)
col.pal = brewer.pal(9,"Blues")

snpCorPlot = snpCor
colnames(snpCorPlot) = rownames(snpCorPlot) =paste0(rse_gene$RNum,":", rse_gene$BrNum)
pdf("qcPlots/pheatmap_ptsd_genotypes.pdf",h=30,w=30)
pheatmap(snpCorPlot, 
		cluster_rows=TRUE, 
		cluster_cols=TRUE,
		color=col.pal)
dev.off()

#### cut dendrogram
dd = as.dist(1-snpCorPlot)
hc = hclust(dd)
pdf("qcPlots/dendrogram_genotypes.pdf",w=150)
palette(brewer.pal(4, "Dark2"))
myplclust(hc, lab.col = as.numeric(factor(rse_gene$Region)))
dev.off()

####### keep checking
snpCor[!upper.tri(snpCor)] = NA

corLong = as.data.frame.table(snpCor)
colnames(corLong) = c("row_RNum", "col_RNum", "genoCor")
corLong = corLong[!is.na(corLong$genoCor),]

corLong$row_BrNum = rse_gene$BrNum[match(corLong$row_RNum, colnames(rse_gene))]
corLong$col_BrNum = rse_gene$BrNum[match(corLong$col_RNum, colnames(rse_gene))]
corLong$row_Region = rse_gene$Region[match(corLong$row_RNum, colnames(rse_gene))]
corLong$col_Region = rse_gene$Region[match(corLong$col_RNum, colnames(rse_gene))]

sigLong = corLong[corLong$genoCor > 0.6,]

checkLong = sigLong[which(sigLong$row_BrNum != sigLong$col_BrNum),]
rownames(checkLong) = NULL

## how many
length(unique(c(checkLong$row_RNum, checkLong$col_RNum))) # 60
length(unique(c(checkLong$row_BrNum, checkLong$col_BrNum))) #21

## view
checkLong[order(checkLong$row_RNum),]

   # row_RNum col_RNum   genoCor row_BrNum col_BrNum row_Region col_Region
# 1    R16748   R16763 0.9274035    Br5956    Br5992   BasoAmyg      DLPFC
# 2    R16749   R16763 0.9454534    Br5956    Br5992      DLPFC      DLPFC
# 3    R16750   R16763 0.9313706    Br5956    Br5992 MedialAmyg      DLPFC
# 4    R16761   R16767 0.6001192    Br5992    Br5812   BasoAmyg MedialAmyg
# 20   R16763   R16334 0.9317288    Br5992    Br5956      DLPFC       dACC
# 5    R16776   R16787 0.9844157    Br5815    Br5972      DLPFC   BasoAmyg
# 6    R16776   R16788 0.9805406    Br5815    Br5972      DLPFC      DLPFC
# 7    R16776   R16789 0.9843563    Br5815    Br5972      DLPFC MedialAmyg
# 19   R16776   R16330 0.9844281    Br5815    Br5972      DLPFC       dACC
# 40   R16796   R16721 0.9186717    Br5980    Br5930   BasoAmyg      DLPFC
# 49   R16796   R16337 0.9345646    Br5980    Br5930   BasoAmyg       dACC
# 38   R16797   R16720 0.9253908    Br5980    Br5930      DLPFC   BasoAmyg
# 42   R16797   R16722 0.9418912    Br5980    Br5930      DLPFC MedialAmyg
# 41   R16798   R16721 0.9428690    Br5980    Br5930 MedialAmyg      DLPFC
# 50   R16798   R16337 0.8950701    Br5980    Br5930 MedialAmyg       dACC
# 21   R16004   R16418 0.9577799    Br5356    Br5443      DLPFC       dACC
# 8    R16100   R17285 0.9406341    Br5415    Br5569      DLPFC MedialAmyg
# 11   R16100   R16129 0.9378679    Br5415    Br5569      DLPFC      DLPFC
# 18   R16100   R16323 0.9130157    Br5415    Br5569      DLPFC       dACC
# 55   R16100   R16488 0.9226106    Br5415    Br5569      DLPFC   BasoAmyg
# 10   R17026   R17034 0.6064900    Br6193    Br8003   BasoAmyg   BasoAmyg
# 13   R17286   R16309 0.9454835    Br5385    Br5400 MedialAmyg       dACC
# 9    R16057   R16078 0.9508336    Br5303    Br5276      DLPFC      DLPFC
# 12   R16057   R16223 0.9694282    Br5303    Br5276      DLPFC       dACC
# 23   R16057   R16437 0.9922977    Br5303    Br5356      DLPFC MedialAmyg
# 30   R16057   R16439 0.9513022    Br5303    Br5276      DLPFC   BasoAmyg
# 34   R16057   R16440 0.9735528    Br5303    Br5276      DLPFC MedialAmyg
# 16   R16078   R16321 0.9126713    Br5276    Br5303      DLPFC       dACC
# 24   R16078   R16437 0.9581174    Br5276    Br5356      DLPFC MedialAmyg
# 27   R16078   R16438 0.9392738    Br5276    Br5303      DLPFC   BasoAmyg
# 14   R16149   R16309 0.9469555    Br5385    Br5400      DLPFC       dACC
# 17   R16223   R16321 0.9021546    Br5276    Br5303       dACC       dACC
# 25   R16223   R16437 0.9692186    Br5276    Br5356       dACC MedialAmyg
# 28   R16223   R16438 0.9578781    Br5276    Br5303       dACC   BasoAmyg
# 15   R16307   R16309 0.9648432    Br5385    Br5400       dACC       dACC
# 56   R16309   R16489 0.9514392    Br5400    Br5385       dACC   BasoAmyg
# 52   R16310   R16487 0.9704922    Br5415    Br5443       dACC MedialAmyg
# 26   R16321   R16437 0.9098143    Br5303    Br5356       dACC MedialAmyg
# 31   R16321   R16439 0.8911596    Br5303    Br5276       dACC   BasoAmyg
# 35   R16321   R16440 0.8984253    Br5303    Br5276       dACC MedialAmyg
# 39   R16329   R16720 0.9372121    Br5980    Br5930       dACC   BasoAmyg
# 43   R16329   R16722 0.9316502    Br5980    Br5930       dACC MedialAmyg
# 22   R16406   R16418 0.9408206    Br5356    Br5443       dACC       dACC
# 44   R16436   R16258 0.9249210    Br5356    Br5443   BasoAmyg      DLPFC
# 51   R16436   R16486 0.9141085    Br5356    Br5443   BasoAmyg   BasoAmyg
# 29   R16437   R16438 0.9651928    Br5356    Br5303 MedialAmyg   BasoAmyg
# 32   R16437   R16439 0.9585352    Br5356    Br5276 MedialAmyg   BasoAmyg
# 36   R16437   R16440 0.9809738    Br5356    Br5276 MedialAmyg MedialAmyg
# 33   R16438   R16439 0.9474103    Br5303    Br5276   BasoAmyg   BasoAmyg
# 37   R16438   R16440 0.9619829    Br5303    Br5276   BasoAmyg MedialAmyg
# 53   R16503   R16487 0.9373975    Br5415    Br5443   BasoAmyg MedialAmyg
# 54   R16504   R16487 0.9540978    Br5415    Br5443 MedialAmyg MedialAmyg
# 45   R16600   R16271 0.9650453    Br5853    Br5753   BasoAmyg      DLPFC
# 46   R16601   R16271 0.9764151    Br5853    Br5753       dACC      DLPFC
# 47   R16602   R16271 0.9724019    Br5853    Br5753 MedialAmyg      DLPFC
# 48   R16275   R16284 0.8224331    Br5806    Br5544      DLPFC      DLPFC
# 57   R16278   R16547 0.9069865    Br5853    Br5753      DLPFC   BasoAmyg
# 58   R16278   R16548 0.9479165    Br5853    Br5753      DLPFC       dACC
# 59   R16278   R16549 0.9441816    Br5853    Br5753      DLPFC MedialAmyg


##########################################
### read in actual DNA genotype data #####
##########################################
genotyped = readVcf("../Genotypes/PTSD_LIBD_VA_Genotypes_n326_GenotypingBarcode.vcf")
mm = match(rownames(snps), rownames(genotyped))
genotyped = genotyped[mm[!is.na(mm)],]
snps = snps[!is.na(mm),]
vcf = vcf[!is.na(mm),]

## get obs geno
snpsGeno = geno(genotyped)$GT
snpsGeno[snpsGeno == "."] = NA
snpsGeno[snpsGeno == "0/0"] = 0
snpsGeno[snpsGeno == "0/1"] = 1
snpsGeno[snpsGeno == "1/1"] = 2
class(snpsGeno) = "numeric"

## flip
table(rowRanges(vcf)$REF == rowRanges(genotyped)$REF  | 
	rowRanges(vcf)$REF == as.character(unlist(alt(genotyped))) )
toFlip = which(rowRanges(vcf)$REF != rowRanges(genotyped)$REF)
snpsGeno[toFlip,] = 2-snpsGeno[toFlip,] 

## correlate
colnames(snpsGeno)= ss(colnames(snpsGeno), "_")
snpCorObs = cor(snps, snpsGeno, use="pairwise.complete.obs")

## check best
bestCor = data.frame(maxIndex = apply(as.matrix(snpCorObs), 1, which.max),
			maxCor = apply(snpCorObs, 1, max, na.rm=TRUE))
bestCor$RNA_BrNum = rse_gene$BrNum[match(rownames(bestCor), colnames(rse_gene))]
bestCor$DNA_BrNum = colnames(snpCorObs)[bestCor$maxIndex]
bestCor$maxIndex = NULL
bestCor$RNum = rownames(bestCor)
bestCor$RNA_Region = rse_gene$Region[match(bestCor$RNum, colnames(rse_gene))]

bestCor = bestCor[order(bestCor$RNA_BrNum),]
checkCor = bestCor[bestCor$RNA_BrNum != bestCor$DNA_BrNum,]

split(checkCor, checkCor$RNA_BrNum)
split(checkCor, checkCor$DNA_BrNum)


#######################
## drop 14 samples ####
dropRNums = c("R16766", "R16763", "R16309", "R16418", "R16436", 
	"R16437", "R16487", "R16437", "R16487","R16100", "R16776", 
	"R17026", "R16253", "R16275", "R17284", "R16558")

bestCor[bestCor$RNum %in% dropRNums,]
          # maxCor RNA_BrNum DNA_BrNum   RNum RNA_Region
# R17284 0.9709129    Br5303    Br5303 R17284 MedialAmyg
# R16436 0.9520771    Br5356    Br5443 R16436   BasoAmyg
# R16437 0.9648288    Br5356    Br5276 R16437 MedialAmyg
# R16309 0.9821465    Br5400    Br5385 R16309       dACC
# R16100 0.9617218    Br5415    Br5569 R16100      DLPFC
# R16418 0.9693929    Br5443    Br5356 R16418       dACC
# R16487 0.9784638    Br5443    Br5415 R16487 MedialAmyg
# R16558 0.4964426    Br5551    Br5757 R16558   BasoAmyg
# R16253 0.4082759    Br5806    Br5544 R16253   BasoAmyg
# R16275 0.6209154    Br5806    Br5544 R16275      DLPFC
# R16766 0.5344626    Br5812    Br5812 R16766      DLPFC
# R16776 0.9921503    Br5815    Br5972 R16776      DLPFC
# R16763 0.9455658    Br5992    Br5956 R16763      DLPFC
# R17026 0.5722310    Br6193    Br6193 R17026   BasoAmyg
bestCor$DNA_BrNum[bestCor$RNum %in% dropRNums] %in% rse_gene$BrNum

rse_gene_drop = rse_gene[,colnames(rse_gene) %in% dropRNums]
rse_gene = rse_gene[,! colnames(rse_gene) %in% dropRNums]
dim(rse_gene)

as.data.frame(colData(rse_gene_drop))

#######################
## seq metrics ########
#######################
rse_gene$Group = factor(rse_gene$Group, levels = c("Control", "PTSD", "MDD"))
rse_gene$Region[rse_gene$Region == "BasoAmyg"] = "AmygBaso"
rse_gene$Region[rse_gene$Region == "MedialAmyg"] = "AmygMedial"
rse_gene$Region = factor(rse_gene$Region)

table(rse_gene$Group, rse_gene$Region)
g = factor(paste0(rse_gene$Region, ":", rse_gene$Group))

vars = c("ERCCsumLogErr", "overallMapRate","concordMapRate", 
	"mitoRate", "totalAssignedGene","rRNA_rate","RIN")
names(vars) = c("ERCC Spike-Ins", "Overall Map Rate",
	"Corcordant Map Rate", "chrM Map Rate", "Exonic Assignment Rate",
		"rRNA Assignment Rate", "RIN")

pdf("qcPlots/qc_metrics_boxplot_byGroupRegion.pdf", h=7,w=8)
par(mar=c(15,6,2,2), cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4, "Dark2"))
for(i in seq(along=vars)) {
	y = colData(rse_gene)[,vars[i]]
	boxplot(y ~ g, outline = FALSE, ylim = range(y),
		las=3, ylab = names(vars)[i])
	points(y ~ jitter(as.numeric(g),amount=0.15),
		pch = 20+as.numeric(rse_gene$Group),
		bg = rse_gene$Region)
}
dev.off()

pdf("qcPlots/qc_metrics_boxplot_byPlate.pdf", h=5,w=8)
par(mar=c(6,6,2,2), cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4, "Dark2"))
for(i in seq(along=vars)) {
	y = colData(rse_gene)[,vars[i]]
	boxplot(y ~ rse_gene$plate, outline = FALSE, ylim = range(y),
		las=3, ylab = names(vars)[i], xlab= "Plate")
	points(y ~ jitter(rse_gene$plate,amount=0.15),
		pch = 20+as.numeric(rse_gene$Group),
		bg = rse_gene$Region)
}
dev.off()

### add flowcell
pheno = readxl::read_excel("VA_DOD_PhaseI_Final_RNA_Samples_sequencing-info-added.xlsx")
rse_gene$Flowcell_1 = pheno$Flowcell_1[match(colnames(rse_gene), pheno$RNum)]
rse_gene$Flowcell_1 = factor(rse_gene$Flowcell_1,
	levels = unique(pheno$Flowcell_1[order(rse_gene$plate)]))

pdf("qcPlots/qc_metrics_boxplot_byFirstFlowcell.pdf", h=5,w=8)
par(mar=c(10,6,2,2), cex.axis=1.5,cex.lab=1.8)
palette(brewer.pal(4, "Dark2"))
for(i in seq(along=vars)) {
	y = colData(rse_gene)[,vars[i]]
	boxplot(y ~ rse_gene$Flowcell_1, outline = FALSE, ylim = range(y),
		las=3, ylab = names(vars)[i], xlab= "")
	points(y ~ jitter(as.numeric(rse_gene$Flowcell_1),amount=0.15),
		pch = 20+as.numeric(rse_gene$Group),
		bg = rse_gene$Region)
}
dev.off()