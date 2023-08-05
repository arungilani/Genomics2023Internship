library(patchwork)
library(dplyr)
library(Seurat)
library(tidyverse)
library(scCustomize)

###setting up subsets of each ident in Hi-MYC data to individually view each level###

#SHAM set E subset
sham <- subset(set_E, orig.ident == c('SHAM-M505', 'SHAM-M514'))
sham$orig.ident <- droplevels(sham$orig.ident)

#week 2 set E subset
wk2 <- subset(set_E, orig.ident == c('Cx-Wk2-M486', 'Cx-Wk2-M488'))
wk2$orig.ident <- droplevels(wk2$orig.ident)

#week 4 set E subset
wk4 <- subset(set_E, orig.ident == c('Cxwk4-M482', 'Cxwk4-M489'))
wk4$orig.ident <- droplevels(wk4$orig.ident)

#week 7 set E subset
wk7 <- subset(set_E, orig.ident == 'Cxwk7-M502')
wk7$orig.ident <- droplevels(wk7$orig.ident)

#week 10 set E subset
wk10 <- subset(set_E, orig.ident == 'Cxwk10-M509')
wk10$orig.ident <- droplevels(wk10$orig.ident)

#recur set E subset
recur <- subset(set_E, orig.ident == c('Cx-M504-recur', 'Cx-M519-recur'))
recur$orig.ident <- droplevels(recur$orig.ident)

###quality control###

#sham set E
sham[['percent.mt']] <- PercentageFeatureSet(sham, pattern = 'mt.')

head(sham@meta.data, 5)

#week 2 set E
wk2[['percent.mt']] <- PercentageFeatureSet(wk2, pattern = 'mt.')

head(wk2@meta.data, 5)

#week 4 set E
wk4[['percent.mt']] <- PercentageFeatureSet(wk4, pattern = 'mt.')

head(wk4@meta.data, 5)

#week 7 set E
wk7[['percent.mt']] <- PercentageFeatureSet(wk7, pattern = 'mt.')

head(wk7@meta.data, 5)

#week 10 set E
wk10[['percent.mt']] <- PercentageFeatureSet(wk10, pattern = 'mt.')

head(wk10@meta.data, 5)

#recurrence set E
recur[['percent.mt']] <- PercentageFeatureSet(recur, pattern = 'mt.')

head(recur@meta.data, 5)

#violin plot for count and feature QC SHAM set E
sham_feature <- QC_Plots_Feature(seurat_object = sham, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#violin plot for count and feature QC week 2 set E
wk2_feature <- QC_Plots_Feature(seurat_object = wk2, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#violin plot for count and feature QC week 4 set E
wk4_feature <- QC_Plots_Feature(seurat_object = wk4, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#violin plot for count and feature QC week 7 set E
wk7_feature <- QC_Plots_Feature(seurat_object = wk7, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#violin plot for count and feature QC week 10 set E
wk10_feature <- QC_Plots_Feature(seurat_object = wk10, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#violin plot for count and feature QC recur set E
recur_feature <- QC_Plots_Feature(seurat_object = recur, feature = c('nFeature_RNA', 'nCount_RNA'), plot_title = 'nCount_RNA', pt.size = 0,) &
  stat_summary(fun.y = median, geom = 'point', size = 15, colour = 'white', shape = 95) &
  ylim(0,2500)

#combining all the plots into one plot
CombinePlots(plots = list(sham_feature, wk2_feature, wk4_feature, wk7_feature, wk10_feature, recur_feature), legend = 'none')


#violin plot for mitochondrial percent SHAM set E
mito_sham <- QC_Plots_Mito(seurat_object = sham, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#violin plot for mitochondrial percent week 2 set E
mito_wk2 <- QC_Plots_Mito(seurat_object = wk2, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#violin plot for mitochondrial percent week 4 set E
mito_wk4 <- QC_Plots_Mito(seurat_object = wk4, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#violin plot for mitochondrial percent week 7 set E
mito_wk7 <- QC_Plots_Mito(seurat_object = wk7, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#violin plot for mitochondrial percent week 10 set E
mito_wk10 <- QC_Plots_Mito(seurat_object = wk10, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#violin plot for mitochondrial percent recur set E
mito_recur <- QC_Plots_Mito(seurat_object = recur, mito_name = 'percent.mt', plot_median = TRUE, pt.size = 0) &
  ylim(0,30)

#combining all the mitochondrial plots into one plot
CombinePlots(plots = list(mito_sham, mito_wk2, mito_wk4, mito_wk7, mito_wk10, mito_recur), legend = 'none')