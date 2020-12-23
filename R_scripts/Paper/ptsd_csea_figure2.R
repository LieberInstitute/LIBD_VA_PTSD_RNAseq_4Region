#reference code: /users/ajaffe/Lieber/Projects/KeriM/csea_functions.R and /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/R_scripts/bipolar_csea_code_cleaned_new.R

## Genes selected for enrichment testing have p < 0.005
library(lattice)
library(jaffelab)
library(SummarizedExperiment)
library(RColorBrewer)

library(gdata)
library(pSI)
library(pSI.data)
data(mouse)
data(human)

fisher.iteration = function (pSIs, candidate.genes,
        background = "data.set",
        thresholds = c(0.05, 0.01, 0.001, 1e-04))  {

        require(IRanges)

      if (background == "human.mouse") {
                load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/CellType_Specificity/kgx_files.rda")
                # just enforce matching by upper case symbol?
                m = match(toupper(rownames(pSIs)), toupper(kgxref.mouse$geneSymbol))
                pSIs = pSIs[!is.na(m),]
                rownames(pSIs) = toupper(rownames(pSIs))

    } else if (background == "mouse.human") {
                load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/CellType_Specificity/kgx_files.rda")
                # just enforce matching by upper case symbol?
                m = match(toupper(rownames(pSIs)), toupper(kgxref.human$geneSymbol))
                pSIs = pSIs[!is.na(m),]
                rownames(pSIs) = toupper(rownames(pSIs))
    } else if(background == "ortholog") {
                ### via https://www.genenames.org/tools/hcop/
                ### human-mouse ortholog as 15 column table
                ortho = read.delim("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/CellType_Specificity/other_data/human_mouse_hcop_fifteen_column.txt",as.is=TRUE)
                ortho = ortho[ortho$human_symbol != "-",]
                m = match(rownames(pSIs), ortho$mouse_symbol)
                pSIs = pSIs[!is.na(m),]
                sym = ortho$human_symbol[m[!is.na(m)]]
                pSIs = pSIs[!duplicated(sym),]
                rownames(pSIs) = sym[!duplicated(sym)]
        } else {
        stop("`background` must be `data.set`, `human.mouse`, `ortholog` or `mouse.human`")
        }

	candidateList = split(candidate.genes,ifelse(toupper(candidate.genes) %in% rownames(pSIs),"Present", "Absent"))

        ## do fisher
        theCandidates = candidateList[["Present"]]
        N = nrow(pSIs)
        names(thresholds) = paste0("psi<",thresholds)
    	statList = lapply(thresholds, function(Th) {
                psiFilter = pSIs < Th # logical matrix based on threshold
                # get symbol names per cell type
                filterList = apply(psiFilter, 2, function(z) names(which(z)))
                inset =  sapply(filterList, function(z) z[z %in% theCandidates])
                # make 2x2 tables
                topleft = sapply(inset, length)
                topright = sapply(filterList, length) - topleft
                bottomleft = length(theCandidates) - topleft
                bottomright = N - topleft - topright - bottomleft

                # the 2x2 table as a matrix
                mat = rbind(topleft,bottomleft,topright,bottomright)
                matList = lapply(as.data.frame(mat), matrix, nrow = 2, nc= 2, byrow=FALSE)
                # fisher p-value
                pv = as.numeric(sapply(matList, function(z) fisher.test(z, alternative = "two.sided")$p.value))
                # and odds ratio for enrichment
                or = sapply(matList, function(z) z[1,1]*z[2,2]/z[1,2]/z[2,1])
                stat = data.frame(OddsRatio = or, pvalue=pv,pvalBH = p.adjust(pv, "BH")) # adjust with BH
                # create IRanges::DataFrame to keep symbols
                stat = DataFrame(stat)
                stat$SymbolsIn = CharacterList(inset)
                return(stat)
        })

	## psi scores and other things to reproduce paper figures
        pSI_out = pSIs[toupper(theCandidates),]
        pvalMat = sapply(statList, function(x) x$pvalue)
        rownames(pvalMat) = colnames(pSIs)
		pvalBHMat = sapply(statList, function(x) x$pvalBH)
		rownames(pvalBHMat) = colnames(pSIs)
        outList = list(pvalMat = pvalMat, pvalBHMat = pvalBHMat, statList = statList,pSI_out = pSI_out,notPresent = candidateList[["Absent"]])
        return(outList)
}

#################################
## read back in merged stats ####
deg_stats = read.csv("../../csvs/merged_deg_stats_allComparisons.csv.gz",
	as.is=TRUE, row.names=1)
colnames(deg_stats) = gsub("DLPFC", "dlPFC", colnames(deg_stats))


###############################
## input gene sets for PTSD ###
###############################

geneSets_cortex_symbols = lapply(with(deg_stats, 
	list(Cortex_either = Cortex_PValue_PTSD < 0.005,
		Cortex_up = Cortex_PValue_PTSD < 0.005 & Cortex_t_PTSD > 0,
		Cortex_down = Cortex_PValue_PTSD < 0.005 & Cortex_t_PTSD < 0,
		dlPFC_either = dlPFC_PValue_PTSD < 0.005,
		dlPFC_up = dlPFC_PValue_PTSD < 0.005 & dlPFC_t_PTSD > 0,
		dlPFC_down = dlPFC_PValue_PTSD < 0.005 & dlPFC_t_PTSD < 0,	
		dACC_either = dACC_PValue_PTSD < 0.005,
		dACC_up = dACC_PValue_PTSD < 0.005 & dACC_t_PTSD > 0,
		dACC_down = dACC_PValue_PTSD < 0.005 & dACC_t_PTSD < 0)),
			function(ii) deg_stats$Symbol[ii])
lengths(geneSets_cortex_symbols)		
			
geneSets_amygdala_symbols = lapply(with(deg_stats, 
	list(Amygdala_either = Amygdala_PValue_PTSD < 0.005,
		Amygdala_up = Amygdala_PValue_PTSD < 0.005 & Amygdala_t_PTSD > 0,
		Amygdala_down = Amygdala_PValue_PTSD < 0.005 & Amygdala_t_PTSD < 0,
		MeA_either = MeA_PValue_PTSD < 0.005,
		MeA_up = MeA_PValue_PTSD < 0.005 & MeA_t_PTSD > 0,
		MeA_down = MeA_PValue_PTSD < 0.005 & MeA_t_PTSD < 0,
		BLA_either = BLA_PValue_PTSD < 0.005,
		BLA_up = BLA_PValue_PTSD < 0.005 & BLA_t_PTSD > 0,
		BLA_down = BLA_PValue_PTSD < 0.005 & BLA_t_PTSD < 0)),
			function(ii) deg_stats$Symbol[ii])
lengths(geneSets_amygdala_symbols)

geneSets = c(geneSets_cortex_symbols, geneSets_amygdala_symbols)

##################
## run csea ######
##################

csea_output_list = parallel::mclapply(geneSets, function(g) {
	cat(".")
	fisher.iteration(pSIs=mouse$psi.out, candidate.genes = g,
			background="mouse.human")
},mc.cores=6)

## get any significance
csea_output_list_stats = lapply(csea_output_list, function(x) {
	for(i in seq(along=x$statList)) {
		x$statList[[i]]$psi = names(x$statList)[i]
		x$statList[[i]]$cellType = rownames(x$statList[[i]])
	}
	do.call("rbind",x$statList)
})
csea_stats = do.call("rbind", csea_output_list_stats)

csea_stats$Comparison = rep(names(csea_output_list_stats), 
	times = sapply(csea_output_list_stats, nrow))
csea_stats = as.data.frame(csea_stats)
csea_stats = csea_stats[,c("cellType", "Comparison", "psi",
			"OddsRatio","pvalue", "pvalBH", "SymbolsIn")]
csea_stats$SymbolsIn = sapply(csea_stats$SymbolsIn, paste,collapse=";")
csea_stats_order = csea_stats[order(csea_stats$pvalue),]
## bonferoni by csea, 630 total tests per psi cut
csea_stats_order$Bonf_sig = csea_stats_order$pvalue*630 < 0.05
write.csv(csea_stats_order, file = "tables/suppTable_SXX_CSEA_output.csv",
	row.names=FALSE)
	
##############################
## fix names for plotting ####
names(csea_output_list) = paste0(gsub("_", "(", names(csea_output_list)), ")")


## get psi p-value
csea_pval_psi05 = sapply(csea_output_list, function(x) x$pvalMat[,"psi<0.05"])
csea_pval_psi01 = sapply(csea_output_list, function(x) x$pvalMat[,"psi<0.01"])
csea_pval_psi001 = sapply(csea_output_list, function(x) x$pvalMat[,"psi<0.001"])
csea_pval_psi0001 = sapply(csea_output_list, function(x) x$pvalMat[,"psi<1e-04"])

### which signif ?
which(csea_pval_psi05 * (nrow(csea_pval_psi05)*ncol(csea_pval_psi05)) < 0.05,
	arr.ind= TRUE)
which(csea_pval_psi01 * (nrow(csea_pval_psi01)*ncol(csea_pval_psi01)) < 0.05,
	arr.ind= TRUE)
which(csea_pval_psi001 * (nrow(csea_pval_psi001)*ncol(csea_pval_psi001)) < 0.05,
	arr.ind= TRUE)
which(csea_pval_psi0001 * (nrow(csea_pval_psi0001)*ncol(csea_pval_psi0001)) < 0.05,
	arr.ind= TRUE)


### make some heatmaps
## plots
theSeq = seq(0,8,by=0.01)

############################
##### up direction #####
############################

mypal_up = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(length(theSeq)))

pdf("figures/figure2/PTSD_csea_heatmaps_up.pdf",h=10,w=6)
csea_pval_psi05_up = csea_pval_psi05[,grep("up", colnames(csea_pval_psi05))]
colnames(csea_pval_psi05_up) = ss(colnames(csea_pval_psi05_up) , "\\(")
print(levelplot(-log10(t(csea_pval_psi05_up)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
		ylab = "", xlab = "", main = "pSI < 0.05"))

csea_pval_psi01_up = csea_pval_psi01[,grep("up", colnames(csea_pval_psi01))]
colnames(csea_pval_psi01_up) = ss(colnames(csea_pval_psi01_up) , "\\(")
print(levelplot(-log10(t(csea_pval_psi01_up)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 0.01"))

csea_pval_psi001_up = csea_pval_psi001[,grep("up", colnames(csea_pval_psi001))]
colnames(csea_pval_psi001_up) = ss(colnames(csea_pval_psi001_up) , "\\(")
print(levelplot(-log10(t(csea_pval_psi001_up)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 0.001"))
		
csea_pval_psi0001_up = csea_pval_psi0001[,grep("up", colnames(csea_pval_psi0001))]
colnames(csea_pval_psi0001_up) = ss(colnames(csea_pval_psi0001_up) , "\\(")
print(levelplot(-log10(t(csea_pval_psi0001_up)), 
	aspect = "fill", at = theSeq,col.regions = mypal_up,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 1e-4"))
		
dev.off()

############################
##### down direction #####
############################

mypal_down = c("white", colorRampPalette(brewer.pal(9,"YlGnBu"))(length(theSeq)))

pdf("figures/figure2/PTSD_csea_heatmaps_down.pdf",h=10,w=6)
csea_pval_psi05_down = csea_pval_psi05[,grep("down", colnames(csea_pval_psi05))]
colnames(csea_pval_psi05_down) = ss(colnames(csea_pval_psi05_down) , "\\(")
print(levelplot(-log10(t(csea_pval_psi05_down)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
		ylab = "", xlab = "", main = "pSI < 0.05"))

csea_pval_psi01_down = csea_pval_psi01[,grep("down", colnames(csea_pval_psi01))]
colnames(csea_pval_psi01_down) = ss(colnames(csea_pval_psi01_down) , "\\(")
print(levelplot(-log10(t(csea_pval_psi01_down)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 0.01"))

csea_pval_psi001_down = csea_pval_psi001[,grep("down", colnames(csea_pval_psi001))]
colnames(csea_pval_psi001_down) = ss(colnames(csea_pval_psi001_down) , "\\(")
print(levelplot(-log10(t(csea_pval_psi001_down)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 0.001"))
		
csea_pval_psi0001_down = csea_pval_psi0001[,grep("down", colnames(csea_pval_psi0001))]
colnames(csea_pval_psi0001_down) = ss(colnames(csea_pval_psi0001_down) , "\\(")
print(levelplot(-log10(t(csea_pval_psi0001_down)), 
	aspect = "fill", at = theSeq,col.regions = mypal_down,	
	scales=list(x=list(rot=90, cex=2), y=list(cex=1.6)),
	par.settings=list(par.main.text=list(cex=2)),
	ylab = "", xlab = "", main = "pSI < 1e-4"))
		
dev.off()
