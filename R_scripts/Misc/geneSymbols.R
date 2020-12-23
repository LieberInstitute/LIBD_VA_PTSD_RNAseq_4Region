load("../rdas/geneStats_DE_qSVA_mergedRegions_threeGroup.rda")

##onlyPTSD
geneSymbols_BasoAmyg_q0.1_onlyPTSD <- subset(dat, BasoAmyg.adj.P.Val_onlyPTSD < 0.1, c(1,13:16))
dim(geneSymbols_BasoAmyg_q0.1_onlyPTSD)
write.csv(geneSymbols_BasoAmyg_q0.1_onlyPTSD, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_q0.1_onlyPTSD.csv")

geneSymbols_MedialAmyg_q0.1_onlyPTSD <- subset(dat, MedialAmyg.adj.P.Val_onlyPTSD < 0.1, c(1,70:73))
dim(geneSymbols_MedialAmyg_q0.1_onlyPTSD)
write.csv(geneSymbols_MedialAmyg_q0.1_onlyPTSD, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_q0.1_onlyPTSD.csv")

geneSymbols_dACC_q0.1_onlyPTSD <- subset(dat, dACC.adj.P.Val_onlyPTSD < 0.1, c(1,32:35))
dim(geneSymbols_dACC_q0.1_onlyPTSD)
write.csv(geneSymbols_dACC_q0.1_onlyPTSD, file = "csvs/geneSymbols/geneSymbols_dACC_q0.1_onlyPTSD.csv")

geneSymbols_DLPFC_q0.1_onlyPTSD <- subset(dat, DLPFC.adj.P.Val_onlyPTSD < 0.1, c(1,51:54))
dim(geneSymbols_DLPFC_q0.1_onlyPTSD)
write.csv(geneSymbols_DLPFC_q0.1_onlyPTSD, file = "csvs/geneSymbols/geneSymbols_DLPFC_q0.1_onlyPTSD.csv")


######################################################################################################################
##PTSD
geneSymbols_BasoAmyg_q0.1_PTSD <- subset(dat, BasoAmyg.adj.P.Val_PTSD < 0.1, c(1,2:5))
dim(geneSymbols_BasoAmyg_q0.1_PTSD)
write.csv(geneSymbols_BasoAmyg_q0.1_PTSD, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_q0.1_PTSD.csv")

geneSymbols_MedialAmyg_q0.1_PTSD <- subset(dat, MedialAmyg.adj.P.Val_PTSD < 0.1, c(1,59:62))
dim(geneSymbols_MedialAmyg_q0.1_PTSD)
write.csv(geneSymbols_MedialAmyg_q0.1_PTSD, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_q0.1_PTSD.csv")

geneSymbols_dACC_q0.1_PTSD <- subset(dat, dACC.adj.P.Val_PTSD < 0.1, c(1,21:24))
dim(geneSymbols_dACC_q0.1_PTSD)
write.csv(geneSymbols_dACC_q0.1_PTSD, file = "csvs/geneSymbols/geneSymbols_dACC_q0.1_PTSD.csv")

geneSymbols_DLPFC_q0.1_PTSD <- subset(dat, DLPFC.adj.P.Val_PTSD < 0.1, c(1,51:54))
dim(geneSymbols_DLPFC_q0.1_PTSD)
write.csv(geneSymbols_DLPFC_q0.1_PTSD, file = "csvs/geneSymbols/geneSymbols_DLPFC_q0.1_PTSD.csv")


######################################################################################################################
##PTSDvsMDD
geneSymbols_BasoAmyg_q0.1_PTSDvsMDD <- subset(dat, BasoAmyg.adj.P.Val_PTSDvsMDD < 0.1, c(1,17:20))
dim(geneSymbols_BasoAmyg_q0.1_PTSDvsMDD)
write.csv(geneSymbols_BasoAmyg_q0.1_PTSDvsMDD, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_q0.1_PTSDvsMDD.csv")

geneSymbols_MedialAmyg_q0.1_PTSDvsMDD <- subset(dat, MedialAmyg.adj.P.Val_PTSDvsMDD < 0.1, c(1,74:77))
dim(geneSymbols_MedialAmyg_q0.1_PTSDvsMDD)
write.csv(geneSymbols_MedialAmyg_q0.1_PTSDvsMDD, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_q0.1_PTSDvsMDD.csv")

geneSymbols_dACC_q0.1_PTSDvsMDD <- subset(dat, dACC.adj.P.Val_PTSDvsMDD < 0.1, c(1,36:39))
dim(geneSymbols_dACC_q0.1_PTSDvsMDD)
write.csv(geneSymbols_dACC_q0.1_PTSDvsMDD, file = "csvs/geneSymbols/geneSymbols_dACC_q0.1_PTSDvsMDD.csv")

geneSymbols_DLPFC_q0.1_PTSDvsMDD <- subset(dat, DLPFC.adj.P.Val_PTSDvsMDD < 0.1, c(1,55:58))
dim(geneSymbols_DLPFC_q0.1_PTSDvsMDD)
write.csv(geneSymbols_DLPFC_q0.1_PTSDvsMDD, file = "csvs/geneSymbols/geneSymbols_DLPFC_q0.1_PTSDvsMDD.csv")


######################################################################################################################
##ANOVA
geneSymbols_BasoAmyg_q0.1_ANOVA <- subset(dat, BasoAmyg.adj.P.Val_ANOVA < 0.1, c(1,10:12))
dim(geneSymbols_BasoAmyg_q0.1_ANOVA)
write.csv(geneSymbols_BasoAmyg_q0.1_ANOVA, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_q0.1_ANOVA.csv")

geneSymbols_MedialAmyg_q0.1_ANOVA <- subset(dat, MedialAmyg.adj.P.Val_ANOVA < 0.1, c(1,67:69))
dim(geneSymbols_MedialAmyg_q0.1_ANOVA)
write.csv(geneSymbols_MedialAmyg_q0.1_ANOVA, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_q0.1_ANOVA.csv")

geneSymbols_dACC_q0.1_ANOVA <- subset(dat, dACC.adj.P.Val_ANOVA < 0.1, c(1,29:31))
dim(geneSymbols_dACC_q0.1_ANOVA)
write.csv(geneSymbols_dACC_q0.1_ANOVA, file = "csvs/geneSymbols/geneSymbols_dACC_q0.1_ANOVA.csv")

geneSymbols_DLPFC_q0.1_ANOVA <- subset(dat, DLPFC.adj.P.Val_ANOVA < 0.1, c(1,48:50))
dim(geneSymbols_DLPFC_q0.1_ANOVA)
write.csv(geneSymbols_DLPFC_q0.1_ANOVA, file = "csvs/geneSymbols/geneSymbols_DLPFC_q0.1_ANOVA.csv")


######################################################################################################################
##MDD
geneSymbols_BasoAmyg_q0.1_MDD <- subset(dat, BasoAmyg.adj.P.Val_MDD < 0.1, c(1,6:9))
dim(geneSymbols_BasoAmyg_q0.1_MDD)
write.csv(geneSymbols_BasoAmyg_q0.1_MDD, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_q0.1_MDD.csv")

geneSymbols_MedialAmyg_q0.1_MDD <- subset(dat, MedialAmyg.adj.P.Val_MDD < 0.1, c(1,63:66))
dim(geneSymbols_MedialAmyg_q0.1_MDD)
write.csv(geneSymbols_MedialAmyg_q0.1_MDD, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_q0.1_MDD.csv")

geneSymbols_dACC_q0.1_MDD <- subset(dat, dACC.adj.P.Val_MDD < 0.1, c(1,25:28))
dim(geneSymbols_dACC_q0.1_MDD)
write.csv(geneSymbols_dACC_q0.1_MDD, file = "csvs/geneSymbols/geneSymbols_dACC_q0.1_MDD.csv")

geneSymbols_DLPFC_q0.1_MDD <- subset(dat, DLPFC.adj.P.Val_MDD < 0.1, c(1,44:47))
dim(geneSymbols_DLPFC_q0.1_MDD)
write.csv(geneSymbols_DLPFC_q0.1_MDD, file = "csvs/geneSymbols/geneSymbols_DLPFC_q0.1_MDD.csv")
