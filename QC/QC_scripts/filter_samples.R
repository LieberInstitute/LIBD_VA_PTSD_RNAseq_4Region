######
library(rtracklayer)
library(SummarizedExperiment)
library(VariantAnnotation)
library(jaffelab)
library(readxl)
library(RColorBrewer)

## add phenotype data
pd = readxl::read_excel("VA_DOD_PhaseI_Final_RNA_Samples_sequencing-info-added.xlsx")
pd = as.data.frame(pd)
rownames(pd) = pd$RNum
colnames(pd) = ss(colnames(pd), "_")
colnames(pd)[c(3,7,8)] = c("RIN", "Sex","Race")

pd$Region = pd$SubBrRegion
pd$Region[is.na(pd$Region)] = pd$BrainRegion[is.na(pd$Region)]
pd$Region = gsub(" ", "", pd$Region)

pd$Region[pd$Region == "Dacc"] = "dACC"
pd$Region[pd$Region == "BasolateralAmyg"] = "BasoAmyg"

####################################
## drop samples based on genotype ##
dropRNums = c("R16057", "R16321", "R16438", "R16437", "R16765",
	"R16771", "R16772", "R16773", "R15816", "R16769", "R16770",
	"R16763", "R16767", "R16253", "R16275", "R16558", "R16776",
	"R16309", "R16100")
length(dropRNums)

## drop 19
pd = pd[ ! pd$RNum %in% dropRNums,]

####################################
### swap subject data ##############
RNums = pd$RNum

## swap pairs
pd$RNum[RNums == "R16720"] = "R16796"
pd$RNum[RNums == "R16796"] = "R16720"

pd$RNum[RNums == "R16722"] = "R16798"
pd$RNum[RNums == "R16798"] = "R16722"

pd$RNum[RNums == "R16278"] = "R16271"
pd$RNum[RNums == "R16271"] = "R16278"

pd$RNum[RNums == "R16418"] = "R16436"
pd$RNum[RNums == "R16436"] = "R16418"

## update brain number
# pd[which(RNums == "R16309"),c(2,5:10)] = pd[which(pd$BrNum == "Br5385")[1],c(2,5:10)] 
# pd[which(RNums == "R16100"),c(2,5:10)] = pd[which(pd$BrNum == "Br5569")[1],c(2,5:10)] 
pd[which(RNums == "R16487"),c(2,5:10)] = pd[which(pd$BrNum == "Br5415")[1],c(2,5:10)] 

### add column
pd$modifiedIDs = FALSE
pd$modifiedIDs[pd$modifiedIDs %in% c("R16798", "R16796", "R16720", "R16722",
	"R16271",  "R16278", "R16436",  "R16418", "R16309", "R16100", "R16487")] = TRUE
	
####################################
## confirm against genotype data ###
####################################

## read in VCF
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.vcf.gz")
colnames(vcf) = ss(ss(colnames(vcf), "/", 9), "_")
vcf = vcf[,pd$RNum] # match up to correct BAM file

# add rs num
snpMap = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")
oo = findOverlaps(vcf, snpMap, type="equal")
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)]
rowData(vcf)$snpRsNum[is.na(rowData(vcf)$snpRsNum)] = rownames(vcf)[is.na(rowData(vcf)$snpRsNum)]
rownames(vcf) = rowData(vcf)$snpRsNum

## filter
vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
          info(vcf)$VDB >0.1,]
		  
## read in obs
genotyped = readVcf("../Genotypes/PTSD_LIBD_VA_Genotypes_n326_GenotypingBarcode.vcf")

mm = match(rownames(vcf), rownames(genotyped))
vcf = vcf[!is.na(mm),]
genotyped = genotyped[mm[!is.na(mm)],]

## RNA-seq SNPs
snpsCalled = geno(vcf)$GT
snpsCalled[snpsCalled == "."] = 0
snpsCalled[snpsCalled == "0/1"] = 1
snpsCalled[snpsCalled == "1/1"] = 2
class(snpsCalled) = "numeric"		  

## DNA genotypes
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
snpCorObs = cor(snpsCalled, snpsGeno, use="pairwise.complete.obs")

## check best
bestCor = data.frame(maxIndex = apply(as.matrix(snpCorObs), 1, which.max),
			maxCor = apply(snpCorObs, 1, max, na.rm=TRUE))
bestCor$RNA_BrNum = pd$BrNum[match(rownames(bestCor), pd$RNum)]
bestCor$DNA_BrNum = colnames(snpCorObs)[bestCor$maxIndex]
bestCor$maxIndex = NULL
bestCor$RNum = rownames(bestCor)
bestCor$RNA_Region = pd$Region[match(bestCor$RNum, pd$RNum)]

checkCor = bestCor[bestCor$RNA_BrNum != bestCor$DNA_BrNum,]
checkCor
# [1] maxCor     RNA_BrNum  DNA_BrNum  RNum       RNA_Region
# <0 rows> (or 0-length row.names)

####################################
### subset counts ##################
####################################

dir.create("count_data")

## load gene counts
load("preprocessed_data/rse_gene_PTSD_VA_LIBD_PostCombined_n1304.Rdata")
rse_gene = rse_gene[,pd$RNum]
colData(rse_gene) = cbind(pd[, c(1:3, 18, 6:12)], colData(rse_gene))
save(rse_gene, file = paste0("count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n",
	ncol(rse_gene), ".Rdata"))

## load exons counts
load("preprocessed_data/rse_exon_PTSD_VA_LIBD_PostCombined_n1304.Rdata")
rse_exon = rse_exon[,pd$RNum]
colData(rse_exon) = cbind(pd[, c(1:3, 18, 6:12)], colData(rse_exon))
save(rse_exon, file = paste0("count_data/rse_exon_PTSD_VA_LIBD_qcAndAnnotated_n",
	ncol(rse_exon), ".Rdata"))

## load junction counts
load("preprocessed_data/rse_jx_PTSD_VA_LIBD_PostCombined_n1304.Rdata")
rse_jx = rse_jx[,pd$RNum]
colData(rse_jx) = cbind(pd[, c(1:3, 18, 6:12)], colData(rse_jx))
save(rse_jx, file = paste0("count_data/rse_jx_PTSD_VA_LIBD_qcAndAnnotated_n",
	ncol(rse_jx), ".Rdata"))

## load tx counts
load("preprocessed_data/rse_tx_PTSD_VA_LIBD_PostCombined_n1304.Rdata")
rse_tx = rse_tx[,pd$RNum]
colData(rse_tx) = cbind(pd[, c(1:3, 18, 6:12)], colData(rse_tx))
save(rse_tx, file = paste0("count_data/rse_tx_PTSD_VA_LIBD_qcAndAnnotated_n",
	ncol(rse_tx), ".Rdata"))

########################
### PLOTS ######

#######################
## seq metrics ########
#######################
rse_gene$Group = factor(rse_gene$Group, levels = c("Control", "PTSD", "MDD"))
rse_gene$Region[rse_gene$Region == "BasoAmyg"] = "AmygB"
rse_gene$Region[rse_gene$Region == "MedialAmyg"] = "AmygM"
rse_gene$Region = factor(rse_gene$Region)

table(rse_gene$Group, rse_gene$Region)
g = factor(paste0(rse_gene$Region, ":", rse_gene$Group))

vars = c("ERCCsumLogErr", "overallMapRate","concordMapRate", 
	"mitoRate", "totalAssignedGene","rRNA_rate","RIN")
names(vars) = c("ERCC Spike-Ins", "Overall Map Rate",
	"Corcordant Map Rate", "chrM Map Rate", "Exonic Assignment Rate",
		"rRNA Assignment Rate", "RIN")
		
## regression
lapply(vars, function(x) {
	summary(lm(colData(rse_gene)[,x] ~ Region + Group, data = colData(rse_gene)))$coef
})
# $`ERCC Spike-Ins`
                # Estimate Std. Error      t value     Pr(>|t|)
# (Intercept) -12.08888668  0.7772364 -15.55368129 4.296370e-50
# RegionAmygM  -0.50769832  0.8985342  -0.56502949 5.721528e-01
# RegiondACC   -2.21572333  0.8985147  -2.46598450 1.379415e-02
# RegionDLPFC  -1.13154171  0.9006214  -1.25640103 2.092000e-01
# GroupPTSD     0.07738875  0.7807616   0.09911957 9.210589e-01
# GroupMDD     -0.87306997  0.7761988  -1.12480203 2.608840e-01

# $`Overall Map Rate`
                # Estimate  Std. Error     t value     Pr(>|t|)
# (Intercept)  0.741205461 0.006996070 105.9459757 0.000000e+00
# RegionAmygM  0.010856925 0.008087898   1.3423667 1.797153e-01
# RegiondACC   0.043566764 0.008087722   5.3867781 8.529668e-08
# RegionDLPFC  0.060614930 0.008106685   7.4771534 1.404730e-13
# GroupPTSD   -0.004598165 0.007027801  -0.6542821 5.130477e-01
# GroupMDD     0.003844483 0.006986731   0.5502550 5.822406e-01

# $`Corcordant Map Rate`
                # Estimate  Std. Error     t value     Pr(>|t|)
# (Intercept)  0.690385056 0.006765537 102.0443898 0.000000e+00
# RegionAmygM  0.012027572 0.007821387   1.5377800 1.243498e-01
# RegiondACC   0.040094323 0.007821217   5.1263536 3.408852e-07
# RegionDLPFC  0.060500401 0.007839555   7.7173257 2.380586e-14
# GroupPTSD   -0.003972000 0.006796222  -0.5844423 5.590258e-01
# GroupMDD     0.003544696 0.006756505   0.5246345 5.999282e-01

# $`chrM Map Rate`
                 # Estimate  Std. Error    t value      Pr(>|t|)
# (Intercept)  0.0334669974 0.001131877 29.5676992 7.475045e-147
# RegionAmygM  0.0016662626 0.001308521  1.2733938  2.031097e-01
# RegiondACC  -0.0038954139 0.001308493 -2.9770239  2.965368e-03
# RegionDLPFC -0.0029818721 0.001311561 -2.2735296  2.315946e-02
# GroupPTSD   -0.0009865628 0.001137011 -0.8676813  3.857316e-01
# GroupMDD    -0.0014225381 0.001130366 -1.2584757  2.084494e-01

# $`Exonic Assignment Rate`
                 # Estimate  Std. Error     t value     Pr(>|t|)
# (Intercept)  0.4131221413 0.006104197 67.67837171 0.000000e+00
# RegionAmygM  0.0001613501 0.007056837  0.02286437 9.817620e-01
# RegiondACC   0.0247981652 0.007056684  3.51413872 4.565805e-04
# RegionDLPFC  0.0556130507 0.007073229  7.86246955 7.950811e-15
# GroupPTSD   -0.0159743968 0.006131883 -2.60513714 9.290291e-03
# GroupMDD    -0.0086793434 0.006096049 -1.42376546 1.547583e-01

# $`rRNA Assignment Rate`
                 # Estimate   Std. Error    t value     Pr(>|t|)
# (Intercept)  1.284329e-04 1.035960e-05 12.3974790 2.057302e-33
# RegionAmygM -7.265605e-06 1.197635e-05 -0.6066628 5.441824e-01
# RegiondACC  -3.485002e-05 1.197609e-05 -2.9099663 3.677233e-03
# RegionDLPFC -5.436428e-05 1.200417e-05 -4.5287830 6.485732e-06
# GroupPTSD    1.441440e-05 1.040659e-05  1.3851232 1.662563e-01
# GroupMDD     1.104193e-05 1.034577e-05  1.0672897 2.860424e-01

# $RIN
               # Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)  6.97647920 0.05346833 130.478720 0.000000e+00
# RegionAmygM  0.08016268 0.06181276   1.296863 1.949122e-01
# RegiondACC   0.59213597 0.06181142   9.579719 4.848410e-21
# RegionDLPFC  0.91814188 0.06195635  14.819174 5.684706e-46
# GroupPTSD   -0.06778389 0.05371084  -1.262015 2.071734e-01
# GroupMDD    -0.07439446 0.05339695  -1.393234 1.637911e-01

lapply(vars, function(x) {
	as.matrix(anova(lm(colData(rse_gene)[,x] ~ 
		Region + Group, data = colData(rse_gene))))
})
# $`ERCC Spike-Ins`
            # Df      Sum Sq  Mean Sq   F value    Pr(>F)
# Region       3    882.1013 294.0338 2.2621642 0.0795594
# Group        2    239.4397 119.7198 0.9210708 0.3983564
# Residuals 1279 166243.1027 129.9790        NA        NA

# $`Overall Map Rate`
            # Df      Sum Sq     Mean Sq    F value       Pr(>F)
# Region       3  0.76444835 0.254816115 24.1964543 3.125971e-15
# Group        2  0.01524563 0.007622815  0.7238361 4.850871e-01
# Residuals 1279 13.46932107 0.010531135         NA           NA

# $`Corcordant Map Rate`
            # Df      Sum Sq     Mean Sq   F value       Pr(>F)
# Region       3  0.71915067 0.239716890 24.340375 2.555423e-15
# Group        2  0.01206876 0.006034380  0.612719 5.420345e-01
# Residuals 1279 12.59626858 0.009848529        NA           NA

# $`chrM Map Rate`
            # Df       Sum Sq      Mean Sq   F value       Pr(>F)
# Region       3 0.0064675157 0.0021558386 7.8207954 3.566421e-05
# Group        2 0.0004579268 0.0002289634 0.8306168 4.360153e-01
# Residuals 1279 0.3525622848 0.0002756546        NA           NA

# $`Exonic Assignment Rate`
            # Df     Sum Sq     Mean Sq   F value       Pr(>F)
# Region       3  0.6676226 0.222540850 27.757844 2.165455e-17
# Group        2  0.0545837 0.027291850  3.404152 3.353640e-02
# Residuals 1279 10.2540293 0.008017224        NA           NA

# $`rRNA Assignment Rate`
            # Df       Sum Sq      Mean Sq  F value       Pr(>F)
# Region       3 6.072168e-07 2.024056e-07 8.765359 9.357228e-06
# Group        2 4.868833e-08 2.434417e-08 1.054246 3.487575e-01
# Residuals 1279 2.953407e-05 2.309153e-08       NA           NA

# $RIN
            # Df     Sum Sq    Mean Sq   F value       Pr(>F)
# Region       3 181.890280 60.6300932 98.566227 2.137224e-57
# Group        2   1.458435  0.7292176  1.185488 3.059327e-01
# Residuals 1279 786.738943  0.6151204        NA           NA

####### 
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
rse_gene$Flowcell_1 = pd$Flowcell[match(colnames(rse_gene), pd$RNum)]
rse_gene$Flowcell_1 = factor(rse_gene$Flowcell_1,
	levels = unique(rse_gene$Flowcell_1[order(rse_gene$plate)]))

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

## rRNA vs chrM
pdf("qcPlots/chrM_vs_rRNA.pdf")
palette(brewer.pal(4, "Dark2"))
par(mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
plot(rse_gene$mitoRate, rse_gene$rRNA_rate,
	pch = 20+as.numeric(rse_gene$Group),
	bg = rse_gene$Region,
	xlab= names(vars)[4], ylab=names(vars)[6])
legend("topright", levels(rse_gene$Region), pch = 15,col=1:4,
	cex=1.5, pt.cex=2)
dev.off()

################
#### do PCA ####
geneRpkm = recount::getRPKM(rse_gene, "Length")
gIndex = rowMeans(geneRpkm) > 0.2
geneExprs = log2(geneRpkm[gIndex,]+1)
pcaGene = prcomp(t(geneExprs))
pcaVars = getPcaVars(pcaGene)

### pc1 vs 2
pdf("qcPlots/PC1vs2.pdf")
par(mar = c(5,6,2,2),cex.axis=1.9,cex.lab=1.9)
palette(brewer.pal(4,"Dark2"))
plot(pcaGene$x[,1], pcaGene$x[,2], 
	pch = 20+as.numeric(rse_gene$Group),
	bg = rse_gene$Region,
	xlab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	ylab = paste0("PC2: ",pcaVars[2], "% Var Expl"))
dev.off()

## correlation
ccPcaMetrics = cor(pcaGene$x[,1:20], 
	as.data.frame(colData(rse_gene)[,vars]))
tPcaMetrics = ccPcaMetrics/sqrt((1-ccPcaMetrics^2)/(ncol(rse_gene)-2))
pvalPcaMetrics = 2*pt(-abs(tPcaMetrics), df = ncol(rse_gene)-1)

## add anova metrics
regVars = c("Region", "Group", "Race", "plate","Flowcell_1", "Position", "Sex", "AgeDeath")
anovaPval = sapply(regVars, function(x) {
	y = colData(rse_gene)[,x]
	apply(pcaGene$x[,1:20], 2, function(p) {
		anova(lm(p ~ y))$`Pr(>F)`[1]
	})
})

library(lattice)
pdf("qcPlots/qc_metrics_vs_PCs.pdf",w=9)
colnames(pvalPcaMetrics) = names(vars)
logPvalMat = -log10(cbind(pvalPcaMetrics, anovaPval))
logPvalMat[logPvalMat>40] = 40
theSeq = seq(0,40,by=0.1)
my.col <- colorRampPalette(c("white","blue"))(length(theSeq))
print(levelplot(logPvalMat, aspect = "fill", 
	at = theSeq,pretty=TRUE,xlab="",ylab="",
	scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
	panel = panel.levelplot.raster, col.regions = my.col))
dev.off()

## position??
rse_gene$PositionNum = factor(rse_gene$Position, 
	levels = paste0(rep(LETTERS[1:8], times=12),
		rep(c("01", "02", "03", "04","05", "06", "07",
			"08", "09", "10", "11","12"), each=8)))
pdf("qcPlots/qc_metrics_PC1vsPosition.pdf",w=20)
par(mar = c(5,6,2,2), cex.axis=1.3, cex.lab=1.8)
boxplot(pcaGene$x[,1] ~ rse_gene$Position, las=3, 
	ylab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	ylim = range(pcaGene$x[,1]), outline=FALSE)	
points(pcaGene$x[,1] ~ jitter(as.numeric(factor(rse_gene$Position)),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)

boxplot(pcaGene$x[,1] ~ rse_gene$PositionNum, las=3, 
	ylab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	ylim = range(pcaGene$x[,1]), outline=FALSE)	
points(pcaGene$x[,1] ~ jitter(as.numeric(rse_gene$PositionNum),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
	
boxplot(rse_gene$ERCCsumLogErr ~ rse_gene$PositionNum, las=3, 
	ylab = names(vars)[1],
	ylim = range(rse_gene$ERCCsumLogErr), outline=FALSE)	
points(rse_gene$ERCCsumLogErr ~ jitter(as.numeric(rse_gene$PositionNum),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
dev.off()

table(rse_gene$Position, rse_gene$Region)

### PC1 vs stuff
pdf("qcPlots/PC1vsExonicRate.pdf")
par(mar = c(5,6,2,2),cex.axis=1.9,cex.lab=1.9)
palette(brewer.pal(4,"Dark2"))
plot(pcaGene$x[,1], rse_gene$totalAssignedGene, 
	pch = 20+as.numeric(rse_gene$Group),
	bg = rse_gene$Region,
	xlab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	ylab =names(vars)[5])
dev.off()

pdf("qcPlots/PC1vsRegion.pdf",w=9)
par(mar = c(5,6,2,2),cex.axis=1.9,cex.lab=1.9)
palette(brewer.pal(4,"Dark2"))
boxplot(pcaGene$x[,1] ~ rse_gene$Region, ylim = range(pcaGene$x[,1]),
	outline=FALSE, ylab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	xlab = "Region")
points(pcaGene$x[,1] ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
dev.off()

pdf("qcPlots/PC1vsPlate.pdf",w=9)
par(mar = c(5,6,2,2),cex.axis=1.9,cex.lab=1.9)
palette(brewer.pal(4,"Dark2"))
boxplot(pcaGene$x[,1] ~ rse_gene$plate, ylim = range(pcaGene$x[,1]),las=3,
	outline=FALSE, ylab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	xlab = "Processing Plate")
points(pcaGene$x[,1] ~ jitter(rse_gene$plate,amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
dev.off()

pdf("qcPlots/PC1vsFlowcell.pdf",w=9)
par(mar = c(13,6,2,2),cex.axis=1.9,cex.lab=1.9)
palette(brewer.pal(4,"Dark2"))
boxplot(pcaGene$x[,1] ~ rse_gene$Flowcell_1, ylim = range(pcaGene$x[,1]),las=3,
	outline=FALSE, ylab = paste0("PC1: ",pcaVars[1], "% Var Expl"),
	xlab = "Main/First Flowcell")
points(pcaGene$x[,1] ~ jitter(as.numeric(rse_gene$Flowcell_1),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
dev.off()


########################
## region scoring ######
########################

rIndexes = splitit(rse_gene$Region)

## for each cell type, find the genes uniquely differentially expressed
##	by alternating which cell types are coded as 0 and which are coded as 1
tstatList <- lapply(rIndexes, function(i) {
	x <- rep(0, ncol(geneExprs))
	x[i] <- 1
	return(genefilter::rowttests(geneExprs, factor(x)))
})

## take 25 genes per cell type that distinguish it from all others
##	first filtering on p-value and then ranking by effect size
## 	decreasing=FALSE instead of TRUE b/c rowttests is backwards with directionality
numProbes=25
probeList <- lapply(tstatList, function(x) {
	y <- x[which(x[, "p.value"] < 1e-80), ]
	yUp <- y[order(y[, "dm"], decreasing = FALSE),] # signs are swapped
	rownames(yUp)[1:numProbes]
})

## filter to those probes we selected
trainingProbes <- unique(unlist(probeList))
trainingProbes = trainingProbes[!is.na(trainingProbes)]
mergeMarkerExprs <- geneExprs[trainingProbes, ]

## this is just code required for the houseman algorithm
form <- as.formula(sprintf("y ~ %s - 1", paste(levels(rse_gene$Region),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~rse_gene$Region - 1))
colnames(phenoDF) <- ss(colnames(phenoDF), "Region", 2)

## do calibration, Houseman algorithm
mergeMarkerExprsScale = scale(mergeMarkerExprs)
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsScale, 
	pheno = phenoDF, modelFix = form)$coefEsts

coefEstsPlot = coefEsts
rownames(coefEstsPlot) = rowData(rse_gene[rownames(coefEsts),])$Symbol
theSeq2 = seq(-4,4,by=0.02)
my.col2 <- colorRampPalette(c("blue", "white","red"))(length(theSeq2))
regionMeansScale = sapply(splitit(rse_gene$Region), function(ii) rowMeans(mergeMarkerExprsScale[,ii]))
pdf("qcPlots/heatmap_regionalPredictions.pdf", h=5,w=12)
print(levelplot(coefEstsPlot, aspect = "fill", 
	at = theSeq2,pretty=TRUE,xlab="",ylab="",
	scales=list(x=list(rot=90)),
	panel = panel.levelplot.raster, col.regions = my.col2))
dev.off()


## get predictions	
predProps = minfi:::projectCellType(scale(geneExprs[rownames(coefEsts),]), coefEsts)
rse_gene$maxPred = colnames(predProps)[apply(predProps, 1, which.max)]
colData(rse_gene) = cbind(colData(rse_gene), predProps)
rse_gene$maxPredProb = apply(predProps,1,max)
rse_gene$labelProb = NA
for(i in 1:ncol(rse_gene)) {
	rse_gene$labelProb[i] = predProps[i,rse_gene$Region[i]]
}

table(rse_gene$maxPred, rse_gene$Region)
table(ss(rse_gene$maxPred, "g"), ss(as.character(rse_gene$Region),"g"))

checkRegions = which(ss(rse_gene$maxPred, "g") != ss(as.character(rse_gene$Region),"g"))
dfCheck= as.data.frame(colData(rse_gene)[checkRegions,
	c("Region", "maxPred","labelProb","maxPredProb")])
dfCheck

## check a few
scn1a = geneExprs[rownames(rse_gene)[rowData(rse_gene)$Symbol == "SCN1A"],]
pdf("qcPlots/SCN1A_regional_example.pdf",w=5)
par(mar=c(9,6,3,2),cex.lab=  1.8, cex.axis=1.8,cex.main=1.8)
## observed
boxplot(scn1a ~ rse_gene$Region, outline=FALSE,las=3,
	ylim = range(scn1a), ylab="log2(RPKM+1)", main="SCN1A")
points(scn1a ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
## recolored
rse_gene$maxPred = factor(rse_gene$maxPred, levels = levels(rse_gene$Region))
boxplot(scn1a ~ rse_gene$Region, outline=FALSE,ylim = range(scn1a),
	las=3,ylab="log2(RPKM+1)", main="SCN1A")
points(scn1a ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
## replotted
boxplot(scn1a ~ rse_gene$maxPred, outline=FALSE,ylim = range(scn1a),
	ylab="log2(RPKM+1)", main="SCN1A",las=3)
points(scn1a ~ jitter(as.numeric(rse_gene$maxPred),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
dev.off()

gabrq = geneExprs[rownames(rse_gene)[rowData(rse_gene)$Symbol == "GABRQ"],]
pdf("qcPlots/GABRQ_regional_example.pdf",w=5)
par(mar=c(9,6,3,2),cex.lab=  1.8, cex.axis=1.8,cex.main=1.8)
## observed
boxplot(gabrq ~ rse_gene$Region, outline=FALSE,las=3,
	ylim = range(gabrq), ylab="log2(RPKM+1)", main="GABRQ")
points(gabrq ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
## recolored
rse_gene$maxPred = factor(rse_gene$maxPred, levels = levels(rse_gene$Region))
boxplot(gabrq ~ rse_gene$Region, outline=FALSE,ylim = range(gabrq),
	las=3,ylab="log2(RPKM+1)", main="GABRQ")
points(gabrq ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
## replotted
boxplot(gabrq ~ rse_gene$maxPred, outline=FALSE,ylim = range(gabrq),
	ylab="log2(RPKM+1)", main="GABRQ",las=3)
points(gabrq ~ jitter(as.numeric(rse_gene$maxPred),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
dev.off()

sytl5 = geneExprs[rownames(rse_gene)[rowData(rse_gene)$Symbol == "SYTL5"],]
pdf("qcPlots/SYTL5_regional_example.pdf",w=5)
par(mar=c(9,6,3,2),cex.lab=  1.8, cex.axis=1.8,cex.main=1.8)
## observed
boxplot(sytl5 ~ rse_gene$Region, outline=FALSE,las=3,
	ylim = range(sytl5), ylab="log2(RPKM+1)", main="SYTL5")
points(sytl5 ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$Region)
## recolored
rse_gene$maxPred = factor(rse_gene$maxPred, levels = levels(rse_gene$Region))
boxplot(sytl5 ~ rse_gene$Region, outline=FALSE,ylim = range(sytl5),
	las=3,ylab="log2(RPKM+1)", main="SYTL5")
points(sytl5 ~ jitter(as.numeric(rse_gene$Region),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
## replotted
boxplot(sytl5 ~ rse_gene$maxPred, outline=FALSE,ylim = range(sytl5),
	ylab="log2(RPKM+1)", main="SYTL5",las=3)
points(sytl5 ~ jitter(as.numeric(rse_gene$maxPred),amount=0.15),
	pch = 20+as.numeric(rse_gene$Group),bg = rse_gene$maxPred)
dev.off()

#### hcluster #####
dd = dist(t(geneExprs))
hc = hclust(dd)
save(dd,hc,file="count_data/dendrogram_data_n1285.rda")

## just marker genes
load("count_data/dendrogram_data_n1285.rda")
ddMarker = dist(t(scale(geneExprs[rownames(coefEsts),])))
hcMarker = hclust(ddMarker)

pdf("qcPlots/dendrogram_log2RpkmPlus1.pdf",w=150)
palette(brewer.pal(4, "Dark2"))
myplclust(hc, lab.col = as.numeric(factor(rse_gene$Region)))
legend("topleft", levels(rse_gene$Region), col = 1:4, pch =15,nc=4,cex=3)
myplclust(hcMarker, lab.col = as.numeric(factor(rse_gene$Region)))
legend("topleft", levels(rse_gene$Region), col = 1:4, pch =15,nc=4,cex=3)
dev.off()

