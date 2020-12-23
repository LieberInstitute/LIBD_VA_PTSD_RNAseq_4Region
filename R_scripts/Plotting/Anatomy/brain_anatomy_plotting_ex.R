library(rgl)
library(misc3d)
library(neurobase)
if (!requireNamespace("aal")) {
  devtools::install_github("muschellij2/aal")
} else {
  library(aal)
}
if (!requireNamespace("MNITemplate")) {
  devtools::install_github("jfortin1/MNITemplate")
} else {
  library(MNITemplate)
}
img = aal_image()
template = readMNI("Brain", res = "1mm")
cut <- 4500
dtemp <- dim(template)
# All of the sections you can label
labs = aal_get_labels()
# Pick the region of the brain you would like to highlight
cingulate = c(labs$index[grep("Cingulate_Ant_L", labs$name)], labs$index[grep("Cingulate_Ant_R", labs$name)])
mask = remake_img(vec = img %in% cingulate, img = img)
### this would be the ``activation'' or surface you want to render 
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cut, alpha = 0.1, draw = TRUE)
contour3d(mask, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("red") )
### add text
text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top")
text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")
rglwidget()


#my image
#devtools::install_github("muschellij2/brainR")
library(neurobase)
library(brainR)
library(abind)

setwd("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/rdas/3D_brain")

template <- readNIfTI(system.file("MNI152_T1_1mm_brain.nii.gz",
package = "brainR"), reorient = FALSE)
dtemp <- dim(template)

img_dacc <- readNIfTI("ROI_MNI_V6_1mm.nii", reorient=FALSE)
labs_dacc<-read.csv("AAL3_labels.csv")
labs_dacc_keep<-c(labs_dacc$index[grep("ACC_sup_L", labs_dacc$name)], labs_dacc$index[grep("ACC_sup_R", labs_dacc$name)])
mask1 = remake_img(vec = img_dacc %in% labs_dacc_keep, img = img_dacc)

img_dlpfc <- readNIfTI("HarvardOxford-cort-maxprob-thr25-1mm.nii.gz",reorient=FALSE)
labs_dlpfc<-read.csv("HarvardOxford-Cortical.csv")
labs_dlpfc_keep<-labs_dlpfc$index[grep("Middle Frontal Gyrus", labs_dlpfc$name)]
mask2 = remake_img(vec = img_dlpfc %in% labs_dlpfc_keep, img = img_dlpfc)

img_bla <- readNIfTI("Juelich-maxprob-thr25-1mm.nii.gz",reorient=FALSE)
labs_bla<-read.csv("Juelich.csv")
labs_bla_keep<-c(labs_bla$index[grep("GM Amygdala_laterobasal group R", labs_bla$name)], labs_bla$index[grep("GM Amygdala_laterobasal group L", labs_bla$name)])
mask3 = remake_img(vec = img_bla %in% labs_bla_keep, img = img_bla)

img_mea <- readNIfTI("A424+2mm.nii.gz", reorient=FALSE)
labs_mea<-read.csv("medialamyg_labels.csv")
labs_mea_keep<-c(labs_mea$index[grep("mAmyg_L", labs_mea$name)], labs_mea$index[grep("mAmyg_R", labs_mea$name)])
mask4 = remake_img(vec = img_mea %in% labs_mea_keep, img = img_mea)

save(mask1, mask2, mask3, mask4, template, file="PTSD_regions.rda")


contour3d(template, level = 4500, alpha = 0.5, draw = TRUE)

contour3d(mask1, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("#FFA10D"))
contour3d(mask2, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("#F17EB4"))
contour3d(mask3, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("#3881B0"))
contour3d(mask4, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("#56A255"))
rglwidget()



#my image
#devtools::install_github("muschellij2/brainR")
library(neurobase)
library(brainR)
library(htmlwidgets)
#library(rgl)
#library(misc3d)

setwd("/home/bree/Desktop/Jaffe_Lab/PTSD_cohort/RNAseq_analysis/Paper1_figures/Figure1/fsleyes")

template <- readNIfTI("MNI152_T1_1mm_brain.nii.gz")
dtemp <- dim(template)

mask1 <- readNIfTI("harvardoxford-cortical_prob_Cingulate_Gyrus_anterior_division_roi_resampled.nii.gz")

mask2 <- readNIfTI("harvardoxford-cortical_prob_Middle_Frontal_Gyrus.nii.gz")

mask3 <- readNIfTI("juelich_prob_GM_Amygdala_laterobasal_group_L.nii.gz")

mask4 <- readNIfTI("juelich_prob_GM_Amygdala_laterobasal_group_R.nii.gz")

mask5 <- readNIfTI("juelich_prob_GM_Amygdala_centromedial_group_L_roi_resampled.nii.gz")

mask6 <- readNIfTI("juelich_prob_GM_Amygdala_centromedial_group_R_roi_resampled.nii.gz")

contour3d(template, level = 4500, alpha = 0.35, draw = TRUE)
contour3d(mask1, level = 0.99, alpha = 1, add = TRUE, color=c("#FF34B3"))
contour3d(mask2, level = 0.99, alpha = 1, add = TRUE, color=c("#B23AEE"))
contour3d(mask3, level = 0.99, alpha = 1, add = TRUE, color=c("#228B22"))
contour3d(mask4, level = 0.99, alpha = 1, add = TRUE, color=c("#228B22"))
contour3d(mask5, level = 0.99, alpha = 1, add = TRUE, color=c("#FF7F00"))
contour3d(mask6, level = 0.99, alpha = 1, add = TRUE, color=c("#FF7F00"))
rglwidget()

#library(htmlwidgets)
#saveWidget(rglwidget(),"x.html")

#hmmmmm
brain<-contour3d(template, level = 4500, alpha = 0.4, draw = FALSE)
scene <- list(brain)
## use 0.99 for level of mask - binary
activation1 <- contour3d(mask1, level = 0.99, alpha = 1, add = TRUE, color = c("#66A61E"), draw = FALSE)
activation2 <- contour3d(mask2, level = 0.99, alpha = 1, add = TRUE, color = c("#E7298A"), draw = FALSE)
activation3 <- contour3d(mask3, level = 0.99, alpha = 1, add = TRUE, color = c("#8DA0CB"), draw = FALSE)

## add these triangles to the list
scene <- c(scene, list(activation1), list(activation2), list(activation3))
## make output image names from image names
fnames <- c("brain.stl", gsub(".nii.gz", ".stl", img, fixed = TRUE))
outfile <- "index_4D_stl.html"
## write the html file out with JavaScript checkboxes
write4D(scene = scene, fnames=fnames, outfile = outfile)
browseURL(outfile)




#Code from https://github.com/muschellij2/brainR and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911196/
### Example data courtesy of Daniel Reich
### Each visit is a binary mask of lesions in the brain
template <- readNIfTI(system.file("MNI152_T1_1mm_brain.nii.gz",
package = "brainR"), reorient = FALSE)
brain <- contour3d(template, level = 4500, alpha = 0.8, draw = FALSE)
imgs <- paste("Visit_", 1:5, ".nii.gz", sep = "")
files <- sapply(imgs, system.file, package = "brainR")
scene <- list(brain)
## loop through images and threshold masks
nimgs <- length(imgs) # get number of images
cols <- rainbow(nimgs) # set colors
for (iimg in 1:nimgs) {
mask <- readNIfTI(files[iimg], reorient = FALSE)[,,,1] # read image mask
## use 0.99 for level of mask - binary
activation <- contour3d(mask, level = 0.99, alpha = 1, add = TRUE,
color = cols[iimg], draw = FALSE)
## add these triangles to the list
scene <- c(scene, list(activation))
}
## make output image names from image names
fnames <- c("brain.stl", gsub(".nii.gz", ".stl", imgs, fixed = TRUE))
outfile <- "index_4D_stl.html"
## write the html file out with JavaScript checkboxes
write4D(scene = scene, fnames = fnames, outfile = outfile, standalone = TRUE)
browseURL(outfile)