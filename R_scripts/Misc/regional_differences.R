###
library(edgeR)
library(jaffelab)
library(SummarizedExperiment)

# load counts counts
load("../../../count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")

## and model results
load("../../Data/rdas/General/ModelMatrices/Main/AllRegions/PTSD_qsvs_Regionintxn.Rdata")

## and line up MDS
load("../../Data/rdas/General/PTSD_LIBD_VA_MDSonly_n326.rda")
mds$BrNum = ss(rownames(mds), "_")
mds = mds[match(rse_gene$BrNum,mds$BrNum),-11]

# and grab duplicate cor
load("../../Data/rdas/DEG/AllRegions/duplicateCorrelation/geneLevel_duplicateCorrelation_limma_forDE_PTSD.rda")

#################
## pairwise regional differences in controls
controlIndex=  which(rse_gene$Group == "Control")
rse_gene_control = rse_gene[,controlIndex] 
qSV_mat_control = qSV_mat[controlIndex,]
mod_control = model.matrix(~ 0 + Region + AgeDeath + 
	Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + 
	overallMapRate + ERCCsumLogErr,	data = colData(rse_gene_control))
colnames(mod_control) = gsub("Region", "", colnames(mod_control))
mds_control = mds[controlIndex,1:3]
## combine
mod_control_final = cbind(mod_control, mds_control, qSV_mat_control)

## model
dge_control = DGEList(counts = assays(rse_gene_control)$counts,
	genes = rowData(rse_gene_control))
dge_control = calcNormFactors(dge_control)
vGene_control = voom(dge_control,mod_control_final, plot=FALSE)

## lmer
fit <- lmFit(vGene_control, block = rse_gene_control$BrNum,
        correlation = gene_dupCorr$consensus.correlation)
eb <- eBayes(fit)

## contrasts for regions
region_combs <- combn(colnames(mod_control_final)[1:4], 2)
region_contrasts <- apply(region_combs, 2, function(x) {
    z <- paste(x, collapse = '-')
    makeContrasts(contrasts = z, levels = mod_control_final)
})
rownames(region_contrasts) <- colnames(mod_control_final)
colnames(region_contrasts) <-
    apply(region_combs, 2, paste, collapse = '-')
eb_contrasts <- eBayes(contrasts.fit(fit, region_contrasts))

## Extract the p-values and add the WM comparisons too
pvals_contrasts <- eb_contrasts$p.value
ts_contrasts <- eb_contrasts$t

## Fing the sig ones
data.frame(
    'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
    'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)

                    # FDRsig Pval10.6sig Pval10.8sig
# BasoAmyg-dACC         2467         473         250
# BasoAmyg-DLPFC        2078         264         109
# BasoAmyg-MedialAmyg    912         123          44
# dACC-DLPFC            5244         900         422
# dACC-MedialAmyg       1773         327         157
# DLPFC-MedialAmyg       743         122          42
