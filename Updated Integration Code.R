library(patchwork)
library(dplyr)
library(Seurat)
library(tidyverse)
library(data.table)
library(forcats)


#combining 10 orig.ident levels into 6 levels
Hi_MYC <- subset(set_E)
x <- levels(Hi_MYC)
levels(Hi_MYC$orig.ident) <- fct_collapse(x, sham = c('SHAM-M505', 'SHAM-M514'), wk2 = c('Cx-Wk2-M486', 'Cx-Wk2-M488'), wk4 = c('Cxwk4-M482', 'Cxwk4-M489'), wk7= 'Cxwk7-M502', wk10 = 'Cxwk10-M509', recur = c('Cx-M504-recur', 'Cx-M519-recur'))
levels(Hi_MYC$orig.ident)

#quality control and normalization
Hi_MYC[['percent.mt']] <- PercentageFeatureSet(Hi_MYC, pattern = 'mt.')
head(Hi_MYC@meta.data, 5)

Hi_MYC <- subset(Hi_MYC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#renaming SHAM-M505 and SHAM-M514 so that we can use the splitobject function in the next step
levels(Hi_MYC$orig.ident)[levels(Hi_MYC$orig.ident)=='sham'] <- 'SHAM'
levels(Hi_MYC$orig.ident)
Hi_MYC$orig.ident <- droplevels(Hi_MYC$orig.ident)

#renaming all the other orig.idents so that we can use the splitobject function  below in the next step
levels(Hi_MYC$orig.ident)[levels(Hi_MYC$orig.ident)==c('recur', 'wk2', 'wk10', 'wk4', 'wk7')] <- 'STIM'
levels(Hi_MYC$orig.ident)
Hi_MYC$orig.ident <- droplevels(Hi_MYC$orig.ident)

#splitobject function used to split the dataframe for integration
Hi_MYC.list <- SplitObject(Hi_MYC, split.by = 'orig.ident')

# normalize and identify integration features for each dataset independently
Hi_MYC.list <- lapply(X = Hi_MYC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
})

features <- SelectIntegrationFeatures(object.list = Hi_MYC.list)

#identifying integration anchors and performing integration
immune.anchors <- FindIntegrationAnchors(object.list = Hi_MYC.list, anchor.features = features)

#creates an 'integrated data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

#elbow plot to dermine number of PCAs shoud we use. I chose 9
ElbowPlot(immune.combined)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 9, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:9)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:9)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindAllMarkers(object=immune.combined)
head(nk.markers)

#explore marker genes for each cluster
FeaturePlot(immune.combined, features = c('Apod', 'Mgp', 'Dcn', 'Fbln1', 'Lum', 'Sfrp1'), min.cutoff = 'q6')