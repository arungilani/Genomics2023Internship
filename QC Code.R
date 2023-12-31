library(patchwork)
library(dplyr)
library(Seurat)
library(tidyverse)

#setting up subsets of each ident in Hi-MYC data to individually view each level

set_A <- subset(set_A, orig.ident == c("SHAM-M505", "SHAM-M514"))
levels(set_A$orig.ident)
set_A$orig.ident <- droplevels(set_A$orig.ident)

recur <- subset(set_E, orig.ident == c('Cx-M504-recur', 'Cx-M519-recur'))
recur$orig.ident <- droplevels(recur$orig.ident)

wk2 <- subset(set_E, orig.ident == c('Cx-Wk2-M486', 'Cx-Wk2-M488'))
wk2$orig.ident <- droplevels(wk2$orig.ident)

wk10 <- subset(set_E, orig.ident == 'Cxwk10-M509')
wk10$orig.ident <- droplevels(wk10$orig.ident)

wk4 <- subset(set_E, orig.ident == c('Cxwk4-M482', 'Cxwk4-M489'))
wk4$orig.ident <- droplevels(wk4$orig.ident)

wk7 <- subset(set_E, orig.ident == 'Cxwk7-M502')
wk7$orig.ident <- droplevels(wk7$orig.ident)

sham <- subset(set_E, orig.ident == c('SHAM-M505', 'SHAM-M514'))
sham$orig.ident <- droplevels(sham$orig.ident)


###quality control###

#set A
set_A[['percent.mt']] <- PercentageFeatureSet(set_A, pattern = 'mt.')

head(set_A@meta.data, 5)

#recurrence set E
recur[['percent.mt']] <- PercentageFeatureSet(recur, pattern = 'mt.')

head(recur@meta.data, 5)

#week 2 set E
wk2[['percent.mt']] <- PercentageFeatureSet(wk2, pattern = 'mt.')

head(wk2@meta.data, 5)

#week 10 set E
wk10[['percent.mt']] <- PercentageFeatureSet(wk10, pattern = 'mt.')

head(wk10@meta.data, 5)

#week 4 set E
wk4[['percent.mt']] <- PercentageFeatureSet(wk4, pattern = 'mt.')

head(wk4@meta.data, 5)

#week 7 set E
wk7[['percent.mt']] <- PercentageFeatureSet(wk7, pattern = 'mt.')

head(wk7@meta.data, 5)

#sham set E
sham[['percent.mt']] <- PercentageFeatureSet(sham, pattern = 'mt.')

head(sham@meta.data, 5)

#violin plot QC set A
VlnPlot(set_A, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC recur set E
VlnPlot(recur, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC week 2 set E
VlnPlot(wk2, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC week 10 set E
VlnPlot(wk10, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC week 4 set E
VlnPlot(wk4, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC week 7 set E
VlnPlot(wk7, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)

#violin plot QC SHAM set E
VlnPlot(sham, pt.size = 0, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) &
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95)