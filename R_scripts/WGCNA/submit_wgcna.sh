qsub -V -pe local 8 -l bluejay,mf=6G,h_vmem=10G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_DLPFC.R
qsub -V -pe local 8 -l bluejay,mf=6G,h_vmem=10G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_BasoAmyg.R
qsub -V -pe local 8 -l bluejay,mf=6G,h_vmem=10G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_dACC.R
qsub -V -pe local 8 -l bluejay,mf=6G,h_vmem=10G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_MedialAmyg.R

qsub -V -pe local 8 -l bluejay,mf=8G,h_vmem=12G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_regionintxn_subset_Amyg.R
qsub -V -pe local 8 -l bluejay,mf=8G,h_vmem=12G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_regionintxn_subset_Cortex.R
qsub -V -pe local 8 -l bluejay,mf=10G,h_vmem=14G,h_stack=256M -cwd -b y R CMD BATCH --no-save wgcna_allregions.R