geneStats_BasoAmygall$sig[geneStats_BasoAmygall$P.Value_onlyPTSD < 0.005] <- "Significant"
geneStats_BasoAmygall$sig[geneStats_BasoAmygall$P.Value_onlyPTSD > 0.005] <- "NotSignificant"
geneStats_MedialAmygall$sig[geneStats_MedialAmygall$P.Value_onlyPTSD < 0.005] <- "Significant"
geneStats_MedialAmygall$sig[geneStats_MedialAmygall$P.Value_onlyPTSD > 0.005] <- "NotSignificant"
geneStats_dACCall$sig[geneStats_dACCall$P.Value_onlyPTSD < 0.005] <- "Significant"
geneStats_dACCall$sig[geneStats_dACCall$P.Value_onlyPTSD > 0.005] <- "NotSignificant"
geneStats_DLPFCall$sig[geneStats_DLPFCall$P.Value_onlyPTSD < 0.005] <- "Significant"
geneStats_DLPFCall$sig[geneStats_DLPFCall$P.Value_onlyPTSD > 0.005] <- "NotSignificant"

groupIndexes_BLA=splitit(geneStats_BasoAmygall$sig)
groupIndexes_dACC=splitit(geneStats_dACCall$sig)
groupIndexes_MeA=splitit(geneStats_MedialAmygall$sig)
groupIndexes_DLPFC=splitit(geneStats_DLPFCall$sig)

bla <- sapply(groupIndexes_BLA, function(x) table(factor(geneStats_BasoAmygall$gene_type)[x]))
mea <- sapply(groupIndexes_MeA, function(x) table(factor(geneStats_MedialAmygall$gene_type)[x]))
dacc <- sapply(groupIndexes_dACC, function(x) table(factor(geneStats_dACCall$gene_type)[x]))
dlpfc <- sapply(groupIndexes_DLPFC, function(x) table(factor(geneStats_DLPFCall$gene_type)[x]))

sigbla <- signif(prop.table(bla,1),3)[,2]
sigmea <- signif(prop.table(mea,1),3)[,2]
sigdacc <- signif(prop.table(dacc,1),3)[,2]
sigdlpfc <- signif(prop.table(dlpfc,1),3)[,2]

fin <- cbind(sigbla, sigmea, sigdacc, sigdlpfc)
fin*100