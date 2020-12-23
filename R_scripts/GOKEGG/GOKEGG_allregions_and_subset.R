##Gene ontology and KEGG broken down by broad region

###################################################
#Cortex
###################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Cortex.rda")
geneStatsall <- as.data.frame(geneStatsall)

## overall
p005all_MDD = geneStatsall[,grep("P.Value_MDD", colnames(geneStatsall))] < 0.005
p005all_PTSD = geneStatsall[,grep("P.Value_PTSD\\b", colnames(geneStatsall))] < 0.005
p005all_onlyPTSD = geneStatsall[,grep("P.Value_onlyPTSD", colnames(geneStatsall))] < 0.005
p005all_PTSDvsMDD = geneStatsall[,grep("P.Value_PTSDvsMDD", colnames(geneStatsall))] < 0.005

p005all_rintxnMDD = geneStatsall[,grep("P.Value_rintxnMDD", colnames(geneStatsall))] < 0.005
p005all_rintxnPTSD = geneStatsall[,grep("P.Value_rintxnPTSD", colnames(geneStatsall))] < 0.005
p005all_rintxnonlyPTSD = geneStatsall[,grep("P.Value_rintxnonlyPTSD", colnames(geneStatsall))] < 0.005
p005all_Region = geneStatsall[,grep("P.Value_Region\\b", colnames(geneStatsall))] < 0.005
p005all_Region_onlyPTSD = geneStatsall[,grep("P.Value_Region_onlyPTSD", colnames(geneStatsall))] < 0.005

## by dir
dirIndex_MDD = grepl("P.Value_MDD", colnames(geneStatsall)) 
p005up_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] > 0
p005down_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] < 0

dirIndex_PTSD = grepl("P.Value_PTSD\\b", colnames(geneStatsall)) 
p005up_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] > 0
p005down_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] < 0

dirIndex_onlyPTSD = grepl("P.Value_onlyPTSD", colnames(geneStatsall)) 
p005up_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] > 0
p005down_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] < 0

dirIndex_PTSDvsMDD = grepl("P.Value_PTSDvsMDD", colnames(geneStatsall)) 
p005up_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] > 0
p005down_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] < 0

dirIndex_rintxnMDD = grepl("P.Value_rintxnMDD", colnames(geneStatsall)) 
p005up_rintxnMDD = geneStatsall[,dirIndex_rintxnMDD] < 0.005 & geneStatsall[,grep("logFC_rintxnMDD", colnames(geneStatsall))] > 0
p005down_rintxnMDD = geneStatsall[,dirIndex_rintxnMDD] < 0.005 & geneStatsall[,grep("logFC_rintxnMDD", colnames(geneStatsall))] < 0

dirIndex_rintxnPTSD = grepl("P.Value_rintxnPTSD", colnames(geneStatsall)) 
p005up_rintxnPTSD = geneStatsall[,dirIndex_rintxnPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnPTSD", colnames(geneStatsall))] > 0
p005down_rintxnPTSD = geneStatsall[,dirIndex_rintxnPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnPTSD", colnames(geneStatsall))] < 0

dirIndex_rintxnonlyPTSD = grepl("P.Value_rintxnonlyPTSD", colnames(geneStatsall)) 
p005up_rintxnonlyPTSD = geneStatsall[,dirIndex_rintxnonlyPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnonlyPTSD", colnames(geneStatsall))] > 0
p005down_rintxnonlyPTSD = geneStatsall[,dirIndex_rintxnonlyPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnonlyPTSD", colnames(geneStatsall))] < 0

dirIndex_Region = grepl("P.Value_Region\\b", colnames(geneStatsall)) 
p005up_Region = geneStatsall[,dirIndex_Region] < 0.005 & geneStatsall[,grep("logFC_Region\\b", colnames(geneStatsall))] > 0
p005down_Region = geneStatsall[,dirIndex_Region] < 0.005 & geneStatsall[,grep("logFC_Region\\b", colnames(geneStatsall))] < 0

dirIndex_Region_onlyPTSD = grepl("P.Value_Region_onlyPTSD", colnames(geneStatsall)) 
p005up_Region_onlyPTSD = geneStatsall[,dirIndex_Region_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_Region_onlyPTSD", colnames(geneStatsall))] > 0
p005down_Region_onlyPTSD = geneStatsall[,dirIndex_Region_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_Region_onlyPTSD", colnames(geneStatsall))] < 0

## checks
sum(p005up_MDD + p005down_MDD)
#707
sum(p005all_MDD)
#707

sum(p005up_PTSD + p005down_PTSD)
#577
sum(p005all_PTSD)
#577

sum(p005up_onlyPTSD + p005down_onlyPTSD)
#301
sum(p005all_onlyPTSD)
#301

sum(p005up_PTSDvsMDD + p005down_PTSDvsMDD)
#331
sum(p005all_PTSDvsMDD)
#331

sum(p005up_rintxnMDD + p005down_rintxnMDD)
#309
sum(p005all_rintxnMDD)
#309

sum(p005up_rintxnPTSD + p005down_rintxnPTSD)
#105
sum(p005all_rintxnPTSD)
#105

sum(p005up_rintxnonlyPTSD + p005down_rintxnonlyPTSD)
#79
sum(p005all_rintxnonlyPTSD)
#79

sum(p005up_Region + p005down_Region)
#3538
sum(p005all_Region)
#3538

sum(p005up_Region_onlyPTSD + p005down_Region_onlyPTSD)
#4119
sum(p005all_Region_onlyPTSD)
#4119

## make matrix
geneInclMat = cbind(p005all_MDD, p005up_MDD, p005down_MDD, p005all_PTSD, p005up_PTSD, p005down_PTSD, p005all_onlyPTSD, p005up_onlyPTSD, p005down_onlyPTSD, p005all_PTSDvsMDD, p005up_PTSDvsMDD, 
p005down_PTSDvsMDD,p005all_rintxnMDD, p005up_rintxnMDD, p005down_rintxnMDD, p005all_rintxnPTSD, p005up_rintxnPTSD, p005down_rintxnPTSD, p005all_rintxnonlyPTSD, p005up_rintxnonlyPTSD, 
p005down_rintxnonlyPTSD, p005all_Region, p005up_Region, p005down_Region, p005all_Region_onlyPTSD, p005up_Region_onlyPTSD, p005down_Region_onlyPTSD)

geneList = apply(geneInclMat, 2, function(x) {
        o = geneStatsall$EntrezID[x]
        as.character(o[!is.na(o)])
})

summary(geneList)
#                         Length Class  Mode     
#p005all_MDD               610   -none- character
#p005up_MDD                149   -none- character
#p005down_MDD              461   -none- character
#p005all_PTSD              485   -none- character
#p005up_PTSD               218   -none- character
#p005down_PTSD             267   -none- character
#p005all_onlyPTSD          243   -none- character
#p005up_onlyPTSD           155   -none- character
#p005down_onlyPTSD          88   -none- character
#p005all_PTSDvsMDD         279   -none- character
#p005up_PTSDvsMDD          217   -none- character
#p005down_PTSDvsMDD         62   -none- character
#p005all_rintxnMDD         213   -none- character
#p005up_rintxnMDD          149   -none- character
#p005down_rintxnMDD         64   -none- character
#p005all_rintxnPTSD         53   -none- character
#p005up_rintxnPTSD          35   -none- character
#p005down_rintxnPTSD        18   -none- character
#p005all_rintxnonlyPTSD     34   -none- character
#p005up_rintxnonlyPTSD      18   -none- character
#p005down_rintxnonlyPTSD    16   -none- character
#p005all_Region           3099   -none- character
#p005up_Region            2196   -none- character
#p005down_Region           903   -none- character
#p005all_Region_onlyPTSD  3639   -none- character
#p005up_Region_onlyPTSD   2629   -none- character
#p005down_Region_onlyPTSD 1010   -none- character

# and universe of expressed genes
universe = as.character(geneStatsall$EntrezID[!is.na(geneStatsall$EntrezID)])

length(universe)
#[1] 17219

## run GO and KEGG
go_cortex <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                readable= TRUE)
kegg_cortex <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)

save(go_cortex, kegg_cortex, file="rdas/all_regions/geneSet_threeGroups_qSVA_Cortex.rda")

##Gene ontology and KEGG broken down by broad region

###################################################
#Amyg
###################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_regionintxn_subset_Amyg.rda")
geneStatsall <- as.data.frame(geneStatsall)

## overall
p005all_MDD = geneStatsall[,grep("P.Value_MDD", colnames(geneStatsall))] < 0.005
p005all_PTSD = geneStatsall[,grep("P.Value_PTSD\\b", colnames(geneStatsall))] < 0.005
p005all_onlyPTSD = geneStatsall[,grep("P.Value_onlyPTSD", colnames(geneStatsall))] < 0.005
p005all_PTSDvsMDD = geneStatsall[,grep("P.Value_PTSDvsMDD", colnames(geneStatsall))] < 0.005

p005all_rintxnMDD = geneStatsall[,grep("P.Value_rintxnMDD", colnames(geneStatsall))] < 0.005
p005all_rintxnPTSD = geneStatsall[,grep("P.Value_rintxnPTSD", colnames(geneStatsall))] < 0.005
p005all_rintxnonlyPTSD = geneStatsall[,grep("P.Value_rintxnonlyPTSD", colnames(geneStatsall))] < 0.005
p005all_Region = geneStatsall[,grep("P.Value_Region\\b", colnames(geneStatsall))] < 0.005
p005all_Region_onlyPTSD = geneStatsall[,grep("P.Value_Region_onlyPTSD", colnames(geneStatsall))] < 0.005

## by dir
dirIndex_MDD = grepl("P.Value_MDD", colnames(geneStatsall)) 
p005up_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] > 0
p005down_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] < 0

dirIndex_PTSD = grepl("P.Value_PTSD\\b", colnames(geneStatsall)) 
p005up_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] > 0
p005down_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] < 0

dirIndex_onlyPTSD = grepl("P.Value_onlyPTSD", colnames(geneStatsall)) 
p005up_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] > 0
p005down_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] < 0

dirIndex_PTSDvsMDD = grepl("P.Value_PTSDvsMDD", colnames(geneStatsall)) 
p005up_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] > 0
p005down_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] < 0

dirIndex_rintxnMDD = grepl("P.Value_rintxnMDD", colnames(geneStatsall)) 
p005up_rintxnMDD = geneStatsall[,dirIndex_rintxnMDD] < 0.005 & geneStatsall[,grep("logFC_rintxnMDD", colnames(geneStatsall))] > 0
p005down_rintxnMDD = geneStatsall[,dirIndex_rintxnMDD] < 0.005 & geneStatsall[,grep("logFC_rintxnMDD", colnames(geneStatsall))] < 0

dirIndex_rintxnPTSD = grepl("P.Value_rintxnPTSD", colnames(geneStatsall)) 
p005up_rintxnPTSD = geneStatsall[,dirIndex_rintxnPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnPTSD", colnames(geneStatsall))] > 0
p005down_rintxnPTSD = geneStatsall[,dirIndex_rintxnPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnPTSD", colnames(geneStatsall))] < 0

dirIndex_rintxnonlyPTSD = grepl("P.Value_rintxnonlyPTSD", colnames(geneStatsall)) 
p005up_rintxnonlyPTSD = geneStatsall[,dirIndex_rintxnonlyPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnonlyPTSD", colnames(geneStatsall))] > 0
p005down_rintxnonlyPTSD = geneStatsall[,dirIndex_rintxnonlyPTSD] < 0.005 & geneStatsall[,grep("logFC_rintxnonlyPTSD", colnames(geneStatsall))] < 0

dirIndex_Region = grepl("P.Value_Region\\b", colnames(geneStatsall)) 
p005up_Region = geneStatsall[,dirIndex_Region] < 0.005 & geneStatsall[,grep("logFC_Region\\b", colnames(geneStatsall))] > 0
p005down_Region = geneStatsall[,dirIndex_Region] < 0.005 & geneStatsall[,grep("logFC_Region\\b", colnames(geneStatsall))] < 0

dirIndex_Region_onlyPTSD = grepl("P.Value_Region_onlyPTSD", colnames(geneStatsall)) 
p005up_Region_onlyPTSD = geneStatsall[,dirIndex_Region_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_Region_onlyPTSD", colnames(geneStatsall))] > 0
p005down_Region_onlyPTSD = geneStatsall[,dirIndex_Region_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_Region_onlyPTSD", colnames(geneStatsall))] < 0

## checks
sum(p005up_MDD + p005down_MDD)
#224
sum(p005all_MDD)
#224

sum(p005up_PTSD + p005down_PTSD)
#281
sum(p005all_PTSD)
#281

sum(p005up_onlyPTSD + p005down_onlyPTSD)
#152
sum(p005all_onlyPTSD)
#152

sum(p005up_PTSDvsMDD + p005down_PTSDvsMDD)
#138
sum(p005all_PTSDvsMDD)
#138

sum(p005up_rintxnMDD + p005down_rintxnMDD)
#194
sum(p005all_rintxnMDD)
#194

sum(p005up_rintxnPTSD + p005down_rintxnPTSD)
#162
sum(p005all_rintxnPTSD)
#162

sum(p005up_rintxnonlyPTSD + p005down_rintxnonlyPTSD)
#155
sum(p005all_rintxnonlyPTSD)
#155

sum(p005up_Region + p005down_Region)
#1943
sum(p005all_Region)
#1943

sum(p005up_Region_onlyPTSD + p005down_Region_onlyPTSD)
#2092
sum(p005all_Region_onlyPTSD)
#2092

## make matrix
geneInclMat = cbind(p005all_MDD, p005up_MDD, p005down_MDD, p005all_PTSD, p005up_PTSD, p005down_PTSD, p005all_onlyPTSD, p005up_onlyPTSD, p005down_onlyPTSD, p005all_PTSDvsMDD, p005up_PTSDvsMDD, 
p005down_PTSDvsMDD,p005all_rintxnMDD, p005up_rintxnMDD, p005down_rintxnMDD, p005all_rintxnPTSD, p005up_rintxnPTSD, p005down_rintxnPTSD, p005all_rintxnonlyPTSD, p005up_rintxnonlyPTSD, 
p005down_rintxnonlyPTSD, p005all_Region, p005up_Region, p005down_Region, p005all_Region_onlyPTSD, p005up_Region_onlyPTSD, p005down_Region_onlyPTSD)

geneList = apply(geneInclMat, 2, function(x) {
        o = geneStatsall$EntrezID[x]
        as.character(o[!is.na(o)])
})

summary(geneList)
#                         Length Class  Mode     
#p005all_MDD               176   -none- character
#p005up_MDD                 97   -none- character
#p005down_MDD               79   -none- character
#p005all_PTSD              209   -none- character
#p005up_PTSD                89   -none- character
#p005down_PTSD             120   -none- character
#p005all_onlyPTSD           88   -none- character
#p005up_onlyPTSD            46   -none- character
#p005down_onlyPTSD          42   -none- character
#p005all_PTSDvsMDD          79   -none- character
#p005up_PTSDvsMDD           43   -none- character
#p005down_PTSDvsMDD         36   -none- character
#p005all_rintxnMDD         128   -none- character
#p005up_rintxnMDD           33   -none- character
#p005down_rintxnMDD         95   -none- character
#p005all_rintxnPTSD         96   -none- character
#p005up_rintxnPTSD          39   -none- character
#p005down_rintxnPTSD        57   -none- character
#p005all_rintxnonlyPTSD     96   -none- character
#p005up_rintxnonlyPTSD      50   -none- character
#p005down_rintxnonlyPTSD    46   -none- character
#p005all_Region           1661   -none- character
#p005up_Region             959   -none- character
#p005down_Region           702   -none- character
#p005all_Region_onlyPTSD  1789   -none- character
#p005up_Region_onlyPTSD    941   -none- character
#p005down_Region_onlyPTSD  848   -none- character

# and universe of expressed genes
universe = as.character(geneStatsall$EntrezID[!is.na(geneStatsall$EntrezID)])

length(universe)
#17219

## run GO and KEGG
go_amyg <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                readable= TRUE)
kegg_amyg <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)

save(go_amyg, kegg_amyg, file="rdas/all_regions/geneSet_threeGroups_qSVA_Amyg.rda")

###################################################
#All 4 regions
###################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(jaffelab)
setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq')

## load data
load("rdas/all_regions/geneStats_allcols_DE_qSVA_lowlyexpressedfilter_allregions_threeGroup.rda",verbose=TRUE)
geneStatsall <- as.data.frame(geneStatsall)

## overall
p005all_MDD = geneStatsall[,grep("P.Value_MDD", colnames(geneStatsall))] < 0.005
p005all_PTSD = geneStatsall[,grep("P.Value_PTSD\\b", colnames(geneStatsall))] < 0.005
p005all_onlyPTSD = geneStatsall[,grep("P.Value_onlyPTSD", colnames(geneStatsall))] < 0.005
p005all_PTSDvsMDD = geneStatsall[,grep("P.Value_PTSDvsMDD", colnames(geneStatsall))] < 0.005

## by dir
dirIndex_MDD = grepl("P.Value_MDD", colnames(geneStatsall)) 
p005up_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] > 0
p005down_MDD = geneStatsall[,dirIndex_MDD] < 0.005 & geneStatsall[,grep("logFC_MDD", colnames(geneStatsall))] < 0

dirIndex_PTSD = grepl("P.Value_PTSD\\b", colnames(geneStatsall)) 
p005up_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] > 0
p005down_PTSD = geneStatsall[,dirIndex_PTSD] < 0.005 & geneStatsall[,grep("logFC_PTSD\\b", colnames(geneStatsall))] < 0

dirIndex_onlyPTSD = grepl("P.Value_onlyPTSD", colnames(geneStatsall)) 
p005up_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] > 0
p005down_onlyPTSD = geneStatsall[,dirIndex_onlyPTSD] < 0.005 & geneStatsall[,grep("logFC_onlyPTSD", colnames(geneStatsall))] < 0

dirIndex_PTSDvsMDD = grepl("P.Value_PTSDvsMDD", colnames(geneStatsall)) 
p005up_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] > 0
p005down_PTSDvsMDD = geneStatsall[,dirIndex_PTSDvsMDD] < 0.005 & geneStatsall[,grep("logFC_PTSDvsMDD", colnames(geneStatsall))] < 0

## checks
sum(p005up_MDD + p005down_MDD)
#612
sum(p005all_MDD)
#612

sum(p005up_PTSD + p005down_PTSD)
#619
sum(p005all_PTSD)
#619

sum(p005up_onlyPTSD + p005down_onlyPTSD)
#230
sum(p005all_onlyPTSD)
#230

sum(p005up_PTSDvsMDD + p005down_PTSDvsMDD)
#200
sum(p005all_PTSDvsMDD)
#200

## make matrix
geneInclMat = cbind(p005all_MDD, p005up_MDD, p005down_MDD, p005all_PTSD, p005up_PTSD, p005down_PTSD, p005all_onlyPTSD, p005up_onlyPTSD, p005down_onlyPTSD, 
p005all_PTSDvsMDD, p005up_PTSDvsMDD,p005down_PTSDvsMDD)

geneList = apply(geneInclMat, 2, function(x) {
        o = geneStatsall$EntrezID[x]
        as.character(o[!is.na(o)])
})

summary(geneList)
#                   Length Class  Mode     
#p005all_MDD        526    -none- character
#p005up_MDD         251    -none- character
#p005down_MDD       275    -none- character
#p005all_PTSD       517    -none- character
#p005up_PTSD        145    -none- character
#p005down_PTSD      372    -none- character
#p005all_onlyPTSD   163    -none- character
#p005up_onlyPTSD     58    -none- character
#p005down_onlyPTSD  105    -none- character
#p005all_PTSDvsMDD  137    -none- character
#p005up_PTSDvsMDD    52    -none- character
#p005down_PTSDvsMDD  85    -none- character

# and universe of expressed genes
universe = as.character(geneStatsall$EntrezID[!is.na(geneStatsall$EntrezID)])

length(universe)
#[1] 17219

## run GO and KEGG
go_allR <- compareCluster(geneList, fun = "enrichGO",
                universe = universe, OrgDb = org.Hs.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
                readable= TRUE)
kegg_allR <- compareCluster(geneList, fun = "enrichKEGG",
                universe = universe, organism = "hsa", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1)

save(go_allR, kegg_allR, file="rdas/all_regions/geneSet_threeGroups_qSVA_allregions.rda")


