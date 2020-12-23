#estimating cell fractions in PTSD samples
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)
library(minfi)
#not sure about the minfi...

setwd('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/')

load('/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata',verbose=TRUE)
#coefEsts = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.csv",as.is=TRUE, row.names=1)
load("/users/ajaffe/Lieber/Projects/RNAseq/MatchedCellComp/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda",verbose=TRUE)

## estimating cell types in PTSD cohort
yExprs = log2(getRPKM(rse_gene,"Length")+1) # log2(RPKM+1)
yExprs_Scaled = scale(yExprs[rownames(coefEsts),]) # filter to cell type genes, and scale

## then do projection, optionally scaling to 100%...this was more necessary with single cell data as the testing dataset, and not really needed for bulk
PropEsts= minfi:::projectCellType(yExprs_Scaled, coefEsts)
PropEsts_Scaled = prop.table(PropEsts,1) #prop.table by row

#Add in disorder
#need to check...
#as always, there is a better way to do the following...
dx_cells <- as.data.frame(PropEsts_Scaled)
dx_cells$Group <- colData(rse_gene)$Group[match(rownames(PropEsts_Scaled),colData(rse_gene)$RNum)]

test <-as.data.frame(PropEsts_Scaled)
test$Group <- colData(rse_gene)$Group
identical(dx_cells,test,attrib.as.set=TRUE)


MDD <- dx_cells[dx_cells$Group == "MDD",1:6]
PTSD <- dx_cells[dx_cells$Group == "PTSD",1:6]
Control <- dx_cells[dx_cells$Group == "Control",1:6]

colMeans(MDD)
colMeans(PTSD)
colMeans(Control)

t.test(MDD[,1],Control[,1])
t.test(MDD[,2],Control[,2]) #sig
t.test(MDD[,3],Control[,3])
t.test(MDD[,4],Control[,4]) #sig
t.test(MDD[,5],Control[,5])
t.test(MDD[,6],Control[,6]) #sig
