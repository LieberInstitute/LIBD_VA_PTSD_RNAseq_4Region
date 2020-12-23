####
library(jaffelab)
library(lmerTest)
dir.create("figures/figure3")

## read in data
dat_image = read.csv("../../Data/RNAScope/Img.csv",as.is=TRUE)
dat_roi = read.csv("../../Data/RNAScope/ROI.csv",as.is=TRUE)

## drop negative control images
dat_image = dat_image[!grepl("NegCtrl", dat_image$Var1),]
dat_roi = dat_roi[dat_roi$name %in% dat_image$Var1,]

## fix names for genes
colnames(dat_roi) = gsub("520", "_CRHBP", colnames(dat_roi))
colnames(dat_roi) = gsub("570", "_SST", colnames(dat_roi))
colnames(dat_roi) = gsub("620", "_CORT", colnames(dat_roi))
colnames(dat_roi) = gsub("690", "_GAD2", colnames(dat_roi))

## extract other phenotype info
dat_roi$donor = ss(dat_roi$name, "_")
dat_roi$region = ss(dat_roi$name, "_",2)
dat_roi$region[dat_roi$region == "BLA2"] = "BLA"
dat_roi$image = paste0(ss(dat_roi$name, "X_",1),"X")

## proportions of rois covered by dots
dat_roi$PP_CRHBP = dat_roi$MP_CRHBP / dat_roi$RVolume
dat_roi$PP_SST = dat_roi$MP_SST / dat_roi$RVolume
dat_roi$PP_CORT = dat_roi$MP_CORT / dat_roi$RVolume
dat_roi$PP_GAD2 = dat_roi$MP_GAD2 / dat_roi$RVolume

########
### note: looked for "red" - ie cort and sst, and imaged there
########

g = paste0("MD_", c("CORT", "SST", "CRHBP", "GAD2"))

signif(cor(log2(dat_roi[,g]+1)),3)


gPairs = as.data.frame(combn(g,2),stringsAsFactors=FALSE)

coexpList = vector("list", ncol(gPairs))
names(coexpList) = gsub("MD_", "", apply(gPairs,2,paste,collapse="_"))

pdf("figures/figure3/pairwise_RNAscope.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
for(i in 1:ncol(gPairs)) {
	g1 = log2(dat_roi[,gPairs[1,i]]+1)
	g2 = log2(dat_roi[,gPairs[2,i]]+1)
	plot(g1, g2, pch=21,bg="grey",
		ylim = c(0,11), xlim = c(0,11),
		xlab = gsub("MD_", "", gPairs[1,i]),
		ylab = gsub("MD_", "", gPairs[2,i]))
	
	coexpList[[i]] = summary(lmer(g1~ g2 + 
				donor + region + (1|image), data=dat_roi))$coef
	
}
dev.off()
	
t(sapply(coexpList,function(x) x[2,]))

## tables
table(dat_roi$MD_CORT > 4, dat_roi$MD_CRHBP > 4,
	dnn = c("CORT", "CRHBP"))
	
	
table(dat_roi$MD_CORT > 4, dat_roi$MD_SST > 4,
	dnn = c("CORT", "SST"))
summary(lmer(log2(MD_CORT+1) ~ log2(MD_SST+1) + 
	donor + region + (1|image), data=dat_roi))$coef
	
	
table(dat_roi$MD_CORT > 4, dat_roi$MD_GAD2 > 4,
	dnn = c("CORT", "GAD2"))