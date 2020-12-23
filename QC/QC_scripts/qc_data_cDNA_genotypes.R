####
###
library(rtracklayer)
library(SummarizedExperiment)
library(VariantAnnotation)
library(jaffelab)
library(readxl)

## read in
geno = read.csv("dod rna snpmark.csv", as.is=TRUE)
geno = DataFrame(geno)
geno$Genotypes = strsplit(geno$SNPmark, "\n")
geno$Genotypes = lapply(geno$Genotypes, function(x) {
	y = ss(x, ": ", 2)
	names(y) = ss(x, ": ", 1)
	y
})

## take the subset in common
tt = table(unlist(lapply(geno$Genotypes, names)))
keepSnps = names( tt[tt==max(tt)])

genoMat = sapply(geno$Genotypes, function(x) {
	x[x  %in% c("NoCall", "NULL", "Confilict")] = NA
	x[keepSnps]
})
colnames(genoMat) = geno$RNum
genoDat = as.data.frame(genoMat,stringsAsFactors=FALSE)
for(i in 1:ncol(genoDat)) genoDat[,i] = as.character(genoDat[,i])
genoDat[is.null(genoDat)] = NA # fix the nulls

#############
## checks ###
#############
apply(genoDat, 1 , function(x) table(x, useNA="ifany"))
genoDat[genoDat == NULL] = 
## make dendrogram
genoMatFac = as.data.frame(genoMat)
genoMatNum = genoMatFac
for(i in 1:ncol(genoMatNum)) genoMatNum[,i] = as.numeric(genoMatNum[,i])
colnames(genoMatNum) = paste0(colnames(genoMatNum),":", geno$BrNum[match(colnames(genoMatNum), geno$RNum)])

ccObs = cor(genoMatNum, use="comp")
hcObs = hclust(as.dist(1-ccObs))
myplclust(hcObs)