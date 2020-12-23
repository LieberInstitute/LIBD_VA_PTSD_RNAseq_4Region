

#load up

load('rdas/dACC/geneStats_DE_qSVA_dACC_threeGroup.rda')
load('rdas/DLPFC/geneStats_DE_qSVA_DLPFC_threeGroup.rda')
load('rdas/BasoAmyg/geneStats_DE_qSVA_BasoAmyg_threeGroup.rda')
load('rdas/MedialAmyg/geneStats_DE_qSVA_MedialAmyg_threeGroup.rda')
geneStats_MedialAmyg <- geneStats_medialamyg

#BasoAmyg

geneSymbols_BasoAmyg_onlyPTSD_p005 <- geneStats_BasoAmyg[geneStats_BasoAmyg$P.Value_onlyPTSD < 0.005, c("Symbol","logFC_onlyPTSD")]
write.csv(geneSymbols_BasoAmyg_onlyPTSD_p005, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_onlyPTSD_p005.csv")

geneSymbols_BasoAmyg_PTSD_p005 <- geneStats_BasoAmyg[geneStats_BasoAmyg$P.Value_PTSD < 0.005, c("Symbol","logFC_PTSD")]
write.csv(geneSymbols_BasoAmyg_PTSD_p005, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_PTSD_p005.csv")

geneSymbols_BasoAmyg_MDD_p005 <- geneStats_BasoAmyg[geneStats_BasoAmyg$P.Value_MDD < 0.005, c("Symbol","logFC_MDD")]
write.csv(geneSymbols_BasoAmyg_MDD_p005, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_MDD_p005.csv")

geneSymbols_BasoAmyg_PTSDvsMDD_p005 <- geneStats_BasoAmyg[geneStats_BasoAmyg$P.Value_PTSDvsMDD < 0.005, c("Symbol","logFC_PTSDvsMDD")]
write.csv(geneSymbols_BasoAmyg_PTSDvsMDD_p005, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_PTSDvsMDD_p005.csv")

geneSymbols_BasoAmyg_ANOVA_p005 <- geneStats_BasoAmyg[geneStats_BasoAmyg$P.Value_ANOVA < 0.005, "Symbol"]
write.csv(geneSymbols_BasoAmyg_ANOVA_p005, file = "csvs/geneSymbols/geneSymbols_BasoAmyg_ANOVA_p005.csv")


#MedialAmyg

geneSymbols_MedialAmyg_onlyPTSD_p005 <- geneStats_MedialAmyg[geneStats_MedialAmyg$P.Value_onlyPTSD < 0.005, c("Symbol","logFC_onlyPTSD")]
write.csv(geneSymbols_MedialAmyg_onlyPTSD_p005, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_onlyPTSD_p005.csv")

geneSymbols_MedialAmyg_PTSD_p005 <- geneStats_MedialAmyg[geneStats_MedialAmyg$P.Value_PTSD < 0.005, c("Symbol","logFC_PTSD")]
write.csv(geneSymbols_MedialAmyg_PTSD_p005, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_PTSD_p005.csv")

geneSymbols_MedialAmyg_MDD_p005 <- geneStats_MedialAmyg[geneStats_MedialAmyg$P.Value_MDD < 0.005, c("Symbol","logFC_MDD")]
write.csv(geneSymbols_MedialAmyg_MDD_p005, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_MDD_p005.csv")

geneSymbols_MedialAmyg_PTSDvsMDD_p005 <- geneStats_MedialAmyg[geneStats_MedialAmyg$P.Value_PTSDvsMDD < 0.005, c("Symbol","logFC_PTSDvsMDD")]
write.csv(geneSymbols_MedialAmyg_PTSDvsMDD_p005, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_PTSDvsMDD_p005.csv")

geneSymbols_MedialAmyg_ANOVA_p005 <- geneStats_MedialAmyg[geneStats_MedialAmyg$P.Value_ANOVA < 0.005, "Symbol"]
write.csv(geneSymbols_MedialAmyg_ANOVA_p005, file = "csvs/geneSymbols/geneSymbols_MedialAmyg_ANOVA_p005.csv")



#dACC

geneSymbols_dACC_onlyPTSD_p005 <- geneStats_dACC[geneStats_dACC$P.Value_onlyPTSD < 0.005, c("Symbol","logFC_onlyPTSD")]
write.csv(geneSymbols_dACC_onlyPTSD_p005, file = "csvs/geneSymbols/geneSymbols_dACC_onlyPTSD_p005.csv")

geneSymbols_dACC_PTSD_p005 <- geneStats_dACC[geneStats_dACC$P.Value_PTSD < 0.005, c("Symbol","logFC_PTSD")]
write.csv(geneSymbols_dACC_PTSD_p005, file = "csvs/geneSymbols/geneSymbols_dACC_PTSD_p005.csv")

geneSymbols_dACC_MDD_p005 <- geneStats_dACC[geneStats_dACC$P.Value_MDD < 0.005, c("Symbol","logFC_MDD")]
write.csv(geneSymbols_dACC_MDD_p005, file = "csvs/geneSymbols/geneSymbols_dACC_MDD_p005.csv")

geneSymbols_dACC_PTSDvsMDD_p005 <- geneStats_dACC[geneStats_dACC$P.Value_PTSDvsMDD < 0.005, c("Symbol","logFC_PTSDvsMDD")]
write.csv(geneSymbols_dACC_PTSDvsMDD_p005, file = "csvs/geneSymbols/geneSymbols_dACC_PTSDvsMDD_p005.csv")

geneSymbols_dACC_ANOVA_p005 <- geneStats_dACC[geneStats_dACC$P.Value_ANOVA < 0.005, "Symbol"]
write.csv(geneSymbols_dACC_ANOVA_p005, file = "csvs/geneSymbols/geneSymbols_dACC_ANOVA_p005.csv")




#DLPFC

geneSymbols_DLPFC_onlyPTSD_p005 <- geneStats_DLPFC[geneStats_DLPFC$P.Value_onlyPTSD < 0.005, c("Symbol","logFC_onlyPTSD")]
write.csv(geneSymbols_DLPFC_onlyPTSD_p005, file = "csvs/geneSymbols/geneSymbols_DLPFC_onlyPTSD_p005.csv")

geneSymbols_DLPFC_PTSD_p005 <- geneStats_DLPFC[geneStats_DLPFC$P.Value_PTSD < 0.005, c("Symbol","logFC_PTSD")]
write.csv(geneSymbols_DLPFC_PTSD_p005, file = "csvs/geneSymbols/geneSymbols_DLPFC_PTSD_p005.csv")

geneSymbols_DLPFC_MDD_p005 <- geneStats_DLPFC[geneStats_DLPFC$P.Value_MDD < 0.005, c("Symbol","logFC_MDD")]
write.csv(geneSymbols_DLPFC_MDD_p005, file = "csvs/geneSymbols/geneSymbols_DLPFC_MDD_p005.csv")

geneSymbols_DLPFC_PTSDvsMDD_p005 <- geneStats_DLPFC[geneStats_DLPFC$P.Value_PTSDvsMDD < 0.005, c("Symbol","logFC_PTSDvsMDD")]
write.csv(geneSymbols_DLPFC_PTSDvsMDD_p005, file = "csvs/geneSymbols/geneSymbols_DLPFC_PTSDvsMDD_p005.csv")

geneSymbols_DLPFC_ANOVA_p005 <- geneStats_DLPFC[geneStats_DLPFC$P.Value_ANOVA < 0.005, "Symbol"]
write.csv(geneSymbols_DLPFC_ANOVA_p005, file = "csvs/geneSymbols/geneSymbols_DLPFC_ANOVA_p005.csv")
















