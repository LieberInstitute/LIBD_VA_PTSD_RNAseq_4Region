#reference code: /users/ajaffe/Lieber/Projects/KeriM/csea_functions.R and /dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/CSEA_BKB/R_scripts/bipolar_csea_code_cleaned_new.R

## Genes selected for enrichment testing have p < 0.005

library(jaffelab)
library(SummarizedExperiment)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

##############################################################################################################################################################
###################onlyPTSD
##############################################################################################################################################################

##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.
#### specificity.index()
### Online tool (http://genetics.wustl.edu/jdlab/csea-tool-2/), built from source (http://genetics.wustl.edu/jdlab/psi_package/)

#try using more updated hg reference by modifying `fisher.iteration()` to call hg38 instead of hg19 and to change `read.table` to `read.delim` for hg38 to fix error
#to do this, we must write our own function
#note, in our fisher.iteration.pSI.BKB function, we use default fisher.test alternative = "two-sided", yet in the pSI package, they use "greater"

### download data, once
#system("wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/kgXref.txt.gz")
#kgxref.mouse <- read.table(gzfile("kgXref.txt.gz"), sep = "\t",
#        row.names = NULL, header = FALSE, stringsAsFactors = FALSE,
#        quote = "")
#colnames(kgxref.mouse) = c("kgID", "mRNA", "spID", "spDisplayID",
#        "geneSymbol", "refseq", "protAcc", "description", "rfamAcc", "tRNAName")
#system("rm kgXref.txt.gz")

#system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/kgXref.txt.gz")
#kgxref.human <- read.delim(gzfile("kgXref.txt.gz"), sep = "\t",
#        row.names = NULL, header = FALSE, stringsAsFactors = FALSE,
#        quote = "")
#colnames(kgxref.human) = c("kgID", "mRNA", "spID", "spDisplayID",
#        "geneSymbol", "refseq", "protAcc", "description", "rfamAcc", "tRNAName")
#system("rm kgXref.txt.gz")
#save(kgxref.mouse, kgxref.human, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/kgx_files.rda")

fisher.iteration.BKB = function (pSIs, candidate.genes,
        background = "data.set",
        thresholds = c(0.05, 0.01, 0.001, 1e-04))  {

        require(IRanges)

      if (background == "human.mouse") {
                load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/kgx_files.rda")
                # just enforce matching by upper case symbol?
                m = match(toupper(rownames(pSIs)), toupper(kgxref.mouse$geneSymbol))
                pSIs = pSIs[!is.na(m),]
                rownames(pSIs) = toupper(rownames(pSIs))

    } else if (background == "mouse.human") {
                load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/kgx_files.rda")
                # just enforce matching by upper case symbol?
                m = match(toupper(rownames(pSIs)), toupper(kgxref.human$geneSymbol))
                pSIs = pSIs[!is.na(m),]
                rownames(pSIs) = toupper(rownames(pSIs))
    } else if(background == "ortholog") {
                ### via https://www.genenames.org/tools/hcop/
                ### human-mouse ortholog as 15 column table
                ortho = read.delim("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/other_data/human_mouse_hcop_fifteen_column.txt",as.is=TRUE)
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

				
#BasoAmyg 
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
write.csv(sig_BasoAmyg, 'csvs/geneSymbols/pSI_input/geneSymbols_BasoAmyg_p005_onlyPTSD.csv',row.names=FALSE)

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
write.csv(sig_MedialAmyg, 'csvs/geneSymbols/pSI_input/geneSymbols_MedialAmyg_p005_onlyPTSD.csv',row.names=FALSE)

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]
write.csv(sig_dACC, 'csvs/geneSymbols/pSI_input/geneSymbols_dACC_p005_onlyPTSD.csv',row.names=FALSE)

#DLPFC				
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_onlyPTSD < 0.005, "Symbol"]					
write.csv(sig_DLPFC, 'csvs/geneSymbols/pSI_input/geneSymbols_DLPFC_p005_onlyPTSD.csv',row.names=FALSE)

library(gdata)
library(pSI)
library(pSI.data)
data(mouse)
data(human)

#specificity.index(pSI.in, pSI.in.filter, bts = 50, p_max = 0.1,e_min = 0.3, hist = FALSE, SI = FALSE)
#fisher.iteration(pSIs, candidate.genes, background = "data.set", p.adjust=TRUE)

dougherty_onlyPTSD_BLA <- fisher.iteration(mouse$psi.out, sig_BasoAmyg,background="mouse.human",p.adjust=TRUE)
dougherty_onlyPTSD_MeA <- fisher.iteration(mouse$psi.out, sig_MedialAmyg,background="mouse.human",p.adjust=TRUE)
dougherty_onlyPTSD_dACC <- fisher.iteration(mouse$psi.out, sig_dACC,background="mouse.human",p.adjust=TRUE)
dougherty_onlyPTSD_DLPFC <- fisher.iteration(mouse$psi.out, sig_DLPFC,background="mouse.human",p.adjust=TRUE)

save(dougherty_onlyPTSD_BLA, dougherty_onlyPTSD_MeA, dougherty_onlyPTSD_dACC, dougherty_onlyPTSD_DLPFC, 
	file = "rdas/CellType_Specificity/Dougherty2014/pSI_fourregions_p005_onlyPTSD_Dougherty2014.rda")

library(jaffelab)
library(SummarizedExperiment)
library(pSI.data)
data(mouse)
data(human)

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

#BasoAmyg 
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]

#DLPFC				
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_onlyPTSD < 0.005, "Symbol"]					


dougherty_onlyPTSD_BLA <- fisher.iteration.BKB(mouse$psi.out, sig_BasoAmyg,background="mouse.human")
dougherty_onlyPTSD_MeA <- fisher.iteration.BKB(mouse$psi.out, sig_MedialAmyg,background="mouse.human")
dougherty_onlyPTSD_dACC <- fisher.iteration.BKB(mouse$psi.out, sig_dACC,background="mouse.human")
dougherty_onlyPTSD_DLPFC <- fisher.iteration.BKB(mouse$psi.out, sig_DLPFC,background="mouse.human")

save(dougherty_onlyPTSD_BLA, dougherty_onlyPTSD_MeA, dougherty_onlyPTSD_dACC, dougherty_onlyPTSD_DLPFC, file ="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Dougherty2014/Dougherty2014_mouse_onlyPTSD_fisheriterationBKB.rda")


##############################################################################################################################################################

##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
#### snRNA-seq from DLPFC
#### Li determined cell type via PCA,taking the top 25 principal components for tSNE/clustering analyses. Then used gene specificity score (from SpecScore.R) to assign cell type to cluster. Also compared to known gene markers.
### Genes selected for enrichment have p < 0.005
					
#get data for pSI calculation
#see bipolar_csea_code_cleaned_new.R for the generation of these reference pSIs
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)

#BasoAmyg
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)				

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)				

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)				

#DLPFC
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)				
					
#significant results
li_DLPFC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)
li_dACC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)
li_MedialAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)
li_BasoAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)

write.csv(li_BasoAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_onlyPTSD_Li2018.csv",row.names=FALSE)
write.csv(li_MedialAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_onlyPTSD_Li2018.csv",row.names=FALSE)
write.csv(li_dACC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_Li2018.csv",row.names=FALSE)			
write.csv(li_DLPFC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_onlyPTSD_Li2018.csv",row.names=FALSE)

colnames(li_DLPFC) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(li_dACC) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(li_MedialAmyg) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(li_BasoAmyg) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")

library(tibble)
library(dplyr)

li_DLPFC_005 <- li_DLPFC %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

li_dACC_005 <- li_dACC %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

li_MedialAmyg_005 <- li_MedialAmyg %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

li_BasoAmyg_005 <- li_BasoAmyg %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

##############################################################################################################################################################

##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.
#used only frontal cortex from Lake's data
#### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930 frontal cortex -  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz
#### Cluster based on k-nearest neighbors, then tSNE on first 150 principal components, then clusters annotated manually on basis of known markers
### Genes selected for enrichment have p < 0.005

#doing only using Lake's data for the frontal cortex
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda",verbose=TRUE)

#BasoAmyg
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_BasoAmyg)				

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_MedialAmyg)				

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_dACC)				

#DLPFC
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_onlyPTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_DLPFC)				

#significant results
lake_DLPFC <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_DLPFC)
lake_dACC <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_dACC)
lake_MedialAmyg <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_MedialAmyg)
lake_BasoAmyg <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_BasoAmyg)


colnames(lake_DLPFC) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(lake_dACC) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(lake_MedialAmyg) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")
colnames(lake_BasoAmyg) <- c("0.05_adj", "0.01_adj", "0.001_adj", "0.0001_adj")


library(tibble)
library(dplyr)

lake_DLPFC_005 <- lake_DLPFC %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

lake_dACC_005 <- lake_dACC %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

lake_MedialAmyg_005 <- lake_MedialAmyg %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")

lake_BasoAmyg_005 <- lake_BasoAmyg %>%
							rownames_to_column(var = "feature") %>%
							filter_all(any_vars(. < 0.005)) %>%
							column_to_rownames(var = "feature")


##############################################################################################################################################################
#########MDD
##############################################################################################################################################################

##### Xu et al. 2014. Cell Type-Specific Expression Analysis to Identify Putative Cellular Mechanisms for Neurogenetic Disorders.

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')
				
#BasoAmyg 
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_MDD < 0.005, "Symbol"]
write.csv(sig_BasoAmyg, 'csvs/geneSymbols/pSI_input/geneSymbols_BasoAmyg_p005_MDD.csv')

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_MDD < 0.005, "Symbol"]
write.csv(sig_MedialAmyg, 'csvs/geneSymbols/pSI_input/geneSymbols_MedialAmyg_p005_MDD.csv')

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_MDD < 0.005, "Symbol"]
write.csv(sig_dACC, 'csvs/geneSymbols/pSI_input/geneSymbols_dACC_p005_MDD.csv')

#DLPFC					
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_MDD < 0.005, "Symbol"]					
write.csv(sig_DLPFC, 'csvs/geneSymbols/pSI_input/geneSymbols_DLPFC_p005_MDD.csv')

##############################################################################################################################################################

##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
			
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)

#BasoAmyg
#onlyPTSD
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)				

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)				

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)				

#DLPFC
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)				
					
#results
li_DLPFC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)
li_dACC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)
li_MedialAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)
li_BasoAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)

write.csv(li_BasoAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_MDD_Li2018.csv")
write.csv(li_MedialAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_MDD_Li2018.csv")
write.csv(li_dACC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_Li2018.csv")			
write.csv(li_DLPFC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_MDD_Li2018.csv")

##############################################################################################################################################################

##### Lake et al. 2018. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain.

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Lake_2018/Lake_2018_UMI_Raw_FrontalCortex_pSIs.rda",verbose=TRUE)

#BasoAmyg
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_BasoAmyg)				

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_MedialAmyg)				

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_dACC)				

#DLPFC
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_MDD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_DLPFC)				


#significant results
lake_DLPFC <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_DLPFC)
lake_dACC <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_dACC)
lake_MedialAmyg <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_MedialAmyg)
lake_BasoAmyg <- fisher.iteration.pSI.BKB(pSIs_fctx_normcounts_bulked_Lake2018,sig_BasoAmyg)


write.csv(lake_BasoAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_MDD_Lake2018.csv")
write.csv(lake_MedialAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_MDD_Lake2018.csv")
write.csv(lake_dACC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_Lake2018.csv")			
write.csv(lake_DLPFC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_MDD_Lake2018.csv")

##############################################################################################################################################################
######PTSD
##############################################################################################################################################################

##### Li et al. 2018. Integrative functional genomic analysis of human brain development and neuropsychiatric risks.
					
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/CSEA/rdas/Li_2018/Li_2018_UMI_Raw_Subset_pSIs.rda",verbose=TRUE)


#BasoAmyg
#PTSD
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
sig_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_PTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)				

#MedialAmyg
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
sig_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_PTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)				

#dACC
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
sig_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_PTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)				

#DLPFC
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_PTSD < 0.005, "Symbol"]
fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)				
					
#results
li_DLPFC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_DLPFC)
li_dACC <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_dACC)
li_MedialAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_MedialAmyg)
li_BasoAmyg <- fisher.iteration.pSI.BKB(pSIs_umi_normcounts_bulked_Li2018,sig_BasoAmyg)

write.csv(li_BasoAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_PTSD_Li2018.csv")
write.csv(li_MedialAmyg, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_PTSD_Li2018.csv")
write.csv(li_dACC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_Li2018.csv")			
write.csv(li_DLPFC, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_PTSD_Li2018.csv")

##############################################################################################################################################################






















################################################################################################
#using matt's single nuclei data
################################################################################################

## Genes selected for enrichment have p < 0.005

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

##########################################################################				
#DLPFC (single nuclei from Matt vs our bulk DLPFC or dACC)					
##########################################################################

#from source on cluster using pSI package (http://genetics.wustl.edu/jdlab/psi_package/)
library(gdata)
library(jaffelab)
library(SummarizedExperiment)
library('readxl')
library('devtools')

##DLPFC

#onlyPTSD
onlyPTSD_DLPFC <- read.csv('csvs/geneSymbols/pSI_input/geneSymbols_DLPFC_p005_onlyPTSD.csv')
colnames(onlyPTSD_DLPFC) <- NULL
rownames(onlyPTSD_DLPFC) <- NULL
onlyPTSD_DLPFC <- onlyPTSD_DLPFC[,2]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xdlpfcSamples_hvg_nonzero_filtered_normd_k13or8_MNTDec2019.rda",verbose=TRUE)

mattDLPFC_onlyPTSD_DLPFC_k8 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k8.normd, onlyPTSD_DLPFC,background="data.set",p.adjust=TRUE)
mattDLPFC_onlyPTSD_DLPFC_k13 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k13.normd, onlyPTSD_DLPFC,background="data.set",p.adjust=TRUE)

write.csv(mattDLPFC_onlyPTSD_DLPFC_k8, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_onlyPTSD_MattDLPFC2019_k8.csv")
write.csv(mattDLPFC_onlyPTSD_DLPFC_k13, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_onlyPTSD_MattDLPFC2019_k13.csv")

#PTSD					
load('rdas/DLPFC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_DLPFC_threeGroup.rda')
sig_DLPFC <- geneStats_DLPFCall[geneStats_DLPFCall$P.Value_PTSD < 0.005, "Symbol"]					
write.csv(sig_DLPFC, 'csvs/geneSymbols/pSI_input/geneSymbols_DLPFC_p005_PTSD.csv')

PTSD_DLPFC <- read.csv('csvs/geneSymbols/pSI_input/geneSymbols_DLPFC_p005_PTSD.csv')
colnames(PTSD_DLPFC) <- NULL
rownames(PTSD_DLPFC) <- NULL
PTSD_DLPFC <- PTSD_DLPFC[,2]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xdlpfcSamples_hvg_nonzero_filtered_normd_k13or8_MNTDec2019.rda",verbose=TRUE)

mattDLPFC_PTSD_DLPFC_k8 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k8.normd, PTSD_DLPFC,background="data.set",p.adjust=TRUE)
mattDLPFC_PTSD_DLPFC_k13 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k13.normd, PTSD_DLPFC,background="data.set",p.adjust=TRUE)

write.csv(mattDLPFC_PTSD_DLPFC_k8, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_PTSD_MattDLPFC2019_k8.csv")
write.csv(mattDLPFC_PTSD_DLPFC_k13, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_DLPFC_p005_PTSD_MattDLPFC2019_k13.csv")

##dACC

#onlyPTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
onlyPTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xdlpfcSamples_hvg_nonzero_filtered_normd_k13or8_MNTDec2019.rda",verbose=TRUE)

mattDLPFC_onlyPTSD_dACC_k8 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k8.normd, onlyPTSD_dACC,background="data.set",p.adjust=TRUE)
mattDLPFC_onlyPTSD_dACC_k13 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k13.normd, onlyPTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(mattDLPFC_onlyPTSD_dACC_k8, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_MattDLPFC2019_k8.csv")
write.csv(mattDLPFC_onlyPTSD_dACC_k13, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_MattDLPFC2019_k13.csv")

#PTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
PTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_PTSD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xdlpfcSamples_hvg_nonzero_filtered_normd_k13or8_MNTDec2019.rda",verbose=TRUE)

mattDLPFC_PTSD_dACC_k8 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k8.normd, PTSD_dACC,background="data.set",p.adjust=TRUE)
mattDLPFC_PTSD_dACC_k13 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k13.normd, PTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(mattDLPFC_PTSD_dACC_k8, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_MattDLPFC2019_k8.csv")
write.csv(mattDLPFC_PTSD_dACC_k13, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_MattDLPFC2019_k13.csv")

#MDD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
MDD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_MDD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xdlpfcSamples_hvg_nonzero_filtered_normd_k13or8_MNTDec2019.rda",verbose=TRUE)

mattDLPFC_MDD_dACC_k8 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k8.normd, MDD_dACC,background="data.set",p.adjust=TRUE)
mattDLPFC_MDD_dACC_k13 <- fisher.iteration.pSI.BKB(pSI.hvg.PFC.k13.normd, MDD_dACC,background="data.set",p.adjust=TRUE)

write.csv(mattDLPFC_MDD_dACC_k8, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_MattDLPFC2019_k8.csv")
write.csv(mattDLPFC_MDD_dACC_k13, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_MDD_MattDLPFC2019_k13.csv")

##########################################################################				
#Amygdala (single nuclei from Matt vs our bulk BLA or MeA)					
##########################################################################

##BasoAmyg

#onlyPTSD
onlyPTSD_BasoAmyg <- read.csv('csvs/geneSymbols/pSI_input/geneSymbols_BasoAmyg_p005_onlyPTSD.csv')
colnames(onlyPTSD_BasoAmyg) <- NULL
rownames(onlyPTSD_BasoAmyg) <- NULL
onlyPTSD_BasoAmyg <- onlyPTSD_BasoAmyg[,2]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xAmygSamples_hvg_nonzero_filtered_normd_k11or6_MNTDec2019.rda", verbose=TRUE)		 

mattAmyg_onlyPTSD_BasoAmyg_k6 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k6.normd, onlyPTSD_BasoAmyg,background="data.set",p.adjust=TRUE)
mattAmyg_onlyPTSD_BasoAmyg_k11 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k11.normd, onlyPTSD_BasoAmyg,background="data.set",p.adjust=TRUE)

write.csv(mattAmyg_onlyPTSD_BasoAmyg_k6, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_onlyPTSD_MattAmyg2019_k6.csv")
write.csv(mattAmyg_onlyPTSD_BasoAmyg_k11, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_onlyPTSD_MattAmyg2019_k11.csv")

#PTSD					
load('rdas/BasoAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_BasoAmyg_threeGroup.rda')
PTSD_BasoAmyg <- geneStats_BasoAmygall[geneStats_BasoAmygall$P.Value_PTSD < 0.005, "Symbol"]
					
#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xAmygSamples_hvg_nonzero_filtered_normd_k11or6_MNTDec2019.rda", verbose=TRUE)		 

mattAmyg_PTSD_BasoAmyg_k6 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k6.normd, PTSD_BasoAmyg,background="data.set",p.adjust=TRUE)
mattAmyg_PTSD_BasoAmyg_k11 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k11.normd, PTSD_BasoAmyg,background="data.set",p.adjust=TRUE)

write.csv(mattAmyg_PTSD_BasoAmyg_k6, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_PTSD_MattAmyg2019_k6.csv")
write.csv(mattAmyg_PTSD_BasoAmyg_k11, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_BasoAmyg_p005_PTSD_MattAmyg2019_k11.csv")

##MedialAmyg

#onlyPTSD
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
onlyPTSD_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xAmygSamples_hvg_nonzero_filtered_normd_k11or6_MNTDec2019.rda", verbose=TRUE)		 

mattAmyg_onlyPTSD_MedialAmyg_k6 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k6.normd, onlyPTSD_MedialAmyg,background="data.set",p.adjust=TRUE)
mattAmyg_onlyPTSD_MedialAmyg_k11 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k11.normd, onlyPTSD_MedialAmyg,background="data.set",p.adjust=TRUE)

write.csv(mattAmyg_onlyPTSD_MedialAmyg_k6, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_onlyPTSD_MattAmyg2019_k6.csv")
write.csv(mattAmyg_onlyPTSD_MedialAmyg_k11, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_onlyPTSD_MattAmyg2019_k11.csv")

#PTSD					
load('rdas/MedialAmyg/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_MedialAmyg_threeGroup.rda')
PTSD_MedialAmyg <- geneStats_MedialAmygall[geneStats_MedialAmygall$P.Value_PTSD < 0.005, "Symbol"]
					
#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2xAmygSamples_hvg_nonzero_filtered_normd_k11or6_MNTDec2019.rda", verbose=TRUE)		 

mattAmyg_PTSD_MedialAmyg_k6 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k6.normd, PTSD_MedialAmyg,background="data.set",p.adjust=TRUE)
mattAmyg_PTSD_MedialAmyg_k11 <- fisher.iteration.pSI.BKB(pSI.hvg.Amyg.k11.normd, PTSD_MedialAmyg,background="data.set",p.adjust=TRUE)

write.csv(mattAmyg_PTSD_MedialAmyg_k6, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_PTSD_MattAmyg2019_k6.csv")
write.csv(mattAmyg_PTSD_MedialAmyg_k11, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_MedialAmyg_p005_PTSD_MattAmyg2019_k11.csv")

##########################################################################				
#sACC (single nuclei from Matt vs our bulk dACC)					
##########################################################################

##dACC

#onlyPTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
onlyPTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2x-sACCsamples_hvg_nonzero_filtered_normd_k7or5_MNTJan2020.rda",verbose=TRUE)

mattsACC_onlyPTSD_dACC_k5 <- fisher.iteration.pSI.BKB(pSI.hvg.sACC.k5.normd, onlyPTSD_dACC,background="data.set",p.adjust=TRUE)
mattsACC_onlyPTSD_dACC_k7 <- fisher.iteration.pSI.BKB(pSI.hvg.sACC.k7.normd, onlyPTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(mattsACC_onlyPTSD_dACC_k5, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_MattsACC2019_k5.csv")
write.csv(mattsACC_onlyPTSD_dACC_k7, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_MattsACC2019_k7.csv")

#PTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
PTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_PTSD < 0.005, "Symbol"]

#load specificity indices (narrow and broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/pSIs_10xpilot_2x-sACCsamples_hvg_nonzero_filtered_normd_k7or5_MNTJan2020.rda",verbose=TRUE)

mattsACC_PTSD_dACC_k5 <- fisher.iteration.pSI.BKB(pSI.hvg.sACC.k5.normd, PTSD_dACC,background="data.set",p.adjust=TRUE)
mattsACC_PTSD_dACC_k7 <- fisher.iteration.pSI.BKB(pSI.hvg.sACC.k7.normd, PTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(mattsACC_PTSD_dACC_k5, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_MattsACC2019_k5.csv")
write.csv(mattsACC_PTSD_dACC_k7, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_MattsACC2019_k7.csv")


##########################################################################
#Allen data
##########################################################################

library(scrattch)
library(scrattch.io)
library(scrattch.hicat)	
								 
#narrow and broad cell types using Cell Diversity in the Human Cortex data, trimmed means by cluster of gene expression
setwd("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen")
trimmed_means <- read.csv("trimmed_means.csv")
length(trimmed_means[,1])
dim(trimmed_means)									 
rownames(trimmed_means) <- trimmed_means[,1]
trimmed_means <- trimmed_means[,-1]
#get rid of any genes that are not expressed
trimmed_means <- trimmed_means[rowSums(trimmed_means) > 0,]
save(trimmed_means, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/AllenBrain_CellDiversity_HumanCortex_trimmed_means.rda")

#make broad cell types
#take means of broad cell groups (ie concatenate narrow cell types into broad)
t_trimmed_means <- t(trimmed_means)
t_trimmed_means <- as.data.frame(t_trimmed_means)
t_trimmed_means$cellType <- rownames(t_trimmed_means)

t_trimmed_means$cellType <- sub("\\..*", "", t_trimmed_means$cellType)

t_broad_trimmed_means <- aggregate(.~cellType, t_trimmed_means, mean)
rownames(t_broad_trimmed_means) <- t_broad_trimmed_means$cellType
t_broad_trimmed_means <- t_broad_trimmed_means[,-1]
broad_trimmed_means <- t(t_broad_trimmed_means)
broad_trimmed_means <- as.data.frame(broad_trimmed_means)

save(broad_trimmed_means, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/AllenBrain_CellDiversity_HumanCortex_broad_trimmed_means.rda")

#should be 20296x9 after, 20296x120 before

trimmed_means_nonzero <- trimmed_means
trimmed_means_nonzero[trimmed_means_nonzero == "0"] <- NA
save(trimmed_means_nonzero, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/AllenBrain_CellDiversity_HumanCortex_trimmed_means_nonzero.rda")

broad_trimmed_means_nonzero <- broad_trimmed_means
broad_trimmed_means_nonzero[broad_trimmed_means_nonzero == "0"] <- NA
save(broad_trimmed_means_nonzero, file="/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/AllenBrain_CellDiversity_HumanCortex_broad_trimmed_means_nonzero.rda")

#generated pSIs using pSIs_Allen_Cortex.R script
#pSIs stored at /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/pSIs_AllenBrain_CellDiversity_HumanCortex_broad_trimmed_means_k9.rda


##########################################################################
#Allen data (whole cortex, single nuclei) vs our bulk dACC

##dACC

#onlyPTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
onlyPTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_onlyPTSD < 0.005, "Symbol"]

#load specificity indices (broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/pSIs_AllenBrain_CellDiversity_HumanCortex_broad_trimmed_means_k9.rda", verbose=TRUE)

allenCortex_onlyPTSD_dACC_k9 <- fisher.iteration.pSI.BKB(Allen_Cortex_broad_trimmed_means_k9, onlyPTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(allenCortex_onlyPTSD_dACC_k9, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_onlyPTSD_AllenCortex_broad_trimmed_means_k9.csv")

#PTSD					
load('rdas/dACC/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_dACC_threeGroup.rda')
PTSD_dACC <- geneStats_dACCall[geneStats_dACCall$P.Value_PTSD < 0.005, "Symbol"]

#load specificity indices (broad cell types)
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/CSEA/Allen/pSIs_AllenBrain_CellDiversity_HumanCortex_broad_trimmed_means_k9.rda", verbose=TRUE)

allenCortex_PTSD_dACC_k9 <- fisher.iteration.pSI.BKB(Allen_Cortex_broad_trimmed_means_k9, PTSD_dACC,background="data.set",p.adjust=TRUE)

write.csv(allenCortex_PTSD_dACC_k9, file = "csvs/geneSymbols/pSI_output/geneSymbols_results_dACC_p005_PTSD_AllenCortex_broad_trimmed_means_k9.csv")





