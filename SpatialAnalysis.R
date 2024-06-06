# Load required libraries
library(Seurat)
library(data.table)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dittoSeq)
library(tictoc)
library(ComplexHeatmap)
library(cowplot)

# Load spatial transcriptomics data for anterior and posterior middle turbinate
MT_ant_data <- Load10X_Spatial(
  data.dir = "/home/bzawisza/MiddleTurbinate/Mt_ant", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  slice = "slice1", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)

MT_post_data <- Load10X_Spatial(
  data.dir = "/home/bzawisza/MiddleTurbinate/Mt_post", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  slice = "slice1", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)

# Calculate the percentage of mitochondrial genes
MT_ant_data[["percent.mt"]] <- PercentageFeatureSet(MT_ant_data, pattern = "^mt-")
MT_post_data[["percent.mt"]] <- PercentageFeatureSet(MT_post_data, pattern = "^mt-")

# Violin plots for QC metrics
VlnPlot(MT_ant_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

VlnPlot(MT_post_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Scatter plots for QC metrics
plot1 <- FeatureScatter(MT_ant_data, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(MT_ant_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
plot1 + plot2

# Subset the data based on QC metrics
MT_ant_subset <- subset(MT_ant_data, subset = nFeature_Spatial < 6000 & nCount_Spatial < 35000)
MT_post_subset <- subset(MT_post_data, subset = nFeature_Spatial < 3500 & nCount_Spatial < 25000)

# Print the number of samples filtered out
cat("Filter out", ncol(MT_ant_data) - ncol(MT_ant_subset), 
    "samples from anterior data, with", ncol(MT_ant_subset), "samples left.\n")
cat("Filter out", ncol(MT_post_data) - ncol(MT_post_subset), 
    "samples from posterior data, with", ncol(MT_post_subset), "samples left.\n")

# Load additional data
SEPscL310_counts <- Read10X(data.dir ="~/seurat/New.SEP310/SEPscL310comb")
SEPscL310_seuratobj <- CreateSeuratObject(counts = SEPscL310_counts, project = "SEPscL310")

# Create a list of subsets
inte.list <- list(MT_ant_subset, MT_post_subset)

# Run SCTransform normalization
for (i in seq_along(inte.list)) {
  inte.list[[i]] <- SCTransform(inte.list[[i]], assay = "Spatial", verbose = TRUE)
}

# Select integration features and find integration anchors
features <- SelectIntegrationFeatures(object.list = inte.list)
MT.anchors <- FindIntegrationAnchors(object.list = inte.list, anchor.features = features)

# Integrate the data
MT.combined <- IntegrateData(anchorset = MT.anchors, k.weight = 50)
DefaultAssay(MT.combined) <- "integrated"

# Standard workflow for visualization and clustering
MT.combined <- ScaleData(MT.combined, verbose = FALSE)
MT.combined <- RunPCA(MT.combined, npcs = 30, verbose = FALSE)
MT.combined <- RunUMAP(MT.combined, reduction = "pca", dims = 1:30)
MT.combined <- FindNeighbors(MT.combined, reduction = "pca", dims = 1:30)
MT.combined <- FindClusters(MT.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(MT.combined, reduction = "umap", group.by = "orig.ident", label = TRUE)
p2 <- DimPlot(MT.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("SEP310-MT_Ant-Integrated")
p1
p2

# Split UMAP plot by original identity
DimPlot(MT.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Feature and violin plots
FeaturePlot(MT.combined, features = "COL1A1")
VlnPlot(MT.combined, features = "COL1A2")
DotPlot(MT.combined, features = "COL1A2")

# Save cell cluster counts
md <- as.data.table(MT.combined@meta.data)
cellnumbers <- md[, .N, by = c("orig.ident", "seurat_clusters")]
write.csv(cellnumbers, file = "clusterswithOrigIdent.Integrated.csv")

# Mapping Query Data
anchors <- list()
for (i in seq_along(inte.list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = SEPscL310_seuratobj,
    query = inte.list[[i]],
    normalization.method = "SCT",
    k.filter = NA,
    k.anchor = 10,
    reference.reduction = "pca", 
    dims = 1:30
  )
}

# Map query data to reference
for (i in seq_along(inte.list)) {
  inte.list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = inte.list[[i]],
    reference = SEPscL310_seuratobj, 
    refdata = list(celltype.1 = "seurat_clusters"), 
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

# DimPlots for predicted cell types
p1 <- DimPlot(inte.list[[1]], group.by = "predicted.celltype.1", label = TRUE) + NoLegend()
p2 <- DimPlot(inte.list[[2]], group.by = "predicted.celltype.1", label = TRUE) + NoLegend()
pf1 <- DimPlot(SEPscL310_seuratobj, group.by = "seurat_clusters", label = TRUE) + NoLegend()
p1 + pf1 + p2

# Transfer predicted cell type scores and IDs
predictions <- list()
for (i in seq_along(inte.list)) {
  predictions[[i]] <- TransferData(anchorset = anchors[[i]], refdata = SEPscL310_seuratobj$seurat_clusters, dims = 1:30)
  inte.list[[i]] <- AddMetaData(inte.list[[i]], metadata = predictions[[i]])
}

# Spatial plots for predicted cell types
SpatialDimPlot(inte.list[[1]], features = "predicted.celltype.1")

S1 <- SpatialPlot(inte.list[[1]], label = TRUE, label.size = 3, group.by = "predicted.id", repel = TRUE) + ggtitle("MT-ANT") + NoLegend()
S2 <- SpatialPlot(inte.list[[2]], label = TRUE, label.size = 3, group.by = "predicted.id", repel = TRUE) + ggtitle("MT-Post") + NoLegend()
S1 + S2

s1 <- SpatialDimPlot(inte.list[[1]], label.size = 3, group.by = "predicted.celltype.1")
s2 <- SpatialDimPlot(inte.list[[2]], label.size = 3, group.by = "predicted.celltype.1")
s1 + s2

# Differential expression analysis after integration
MT.markers <- FindAllMarkers(MT.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MT.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(MT.markers, file = "~/MT.integrated.csv")
