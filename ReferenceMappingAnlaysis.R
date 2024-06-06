# Required Libraries
library(Seurat)
library(data.table)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dittoSeq)
library(presto)
library(tictoc)
library(ComplexHeatmap)
library(stringr)

# Define Cell Cycle Genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cycle.genes <- c(s.genes, g2m.genes)

# Helper Functions
create_seurat_object <- function(data_dir, project, min_cells = 0, nFeature_RNA_filter = 0, mt_pattern = "^MT-") {
  counts <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = counts, project = project, min.cells = min_cells)
  seurat_obj@meta.data$orig.ident <- project
  if (nFeature_RNA_filter > 0) {
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > nFeature_RNA_filter)
  }
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  return(seurat_obj)
}

process_seurat_object <- function(seurat_obj, vars_to_regress = c("S.Score", "G2M.Score")) {
  seurat_obj <- SCTransform(seurat_obj, verbose = TRUE, vars.to.regress = vars_to_regress)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = 30)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30, return.model = TRUE)
  return(seurat_obj)
}

# Import and Process Datasets
SEPscL310_seuratobj <- create_seurat_object("~/seurat/New.SEP310/SEPscL310comb", "SEPscL310")
SEPscL310_seuratobj <- process_seurat_object(SEPscL310_seuratobj)

SEPscL036_seuratobj <- create_seurat_object("~/seurat/SEPscL036/filtered_feature_bc_matrix", "SEPscL036", min_cells = 3, nFeature_RNA_filter = 2000)
SEPscL181_seuratobj <- create_seurat_object("~/seurat/SEPscL181/filtered_feature_bc_matrix", "SEPscL181", nFeature_RNA_filter = 3500)

# Rename Identifiers
SEPscL036_seuratobj$orig.ident <- "CNON-CTRL"
SEPscL181_seuratobj$orig.ident <- "CNON-SCZ"

# Create Reference List
reference_list <- list(CNON_CTRL = SEPscL036_seuratobj, CNON_SCZ = SEPscL181_seuratobj)

# Process Reference List
for (i in names(reference_list)) {
  reference_list[[i]] <- process_seurat_object(reference_list[[i]])
}

# Plot Mitochondrial Percentage
V1 <- VlnPlot(reference_list$CNON_CTRL, features = "percent.mt") + scale_y_continuous(limits = c(0, 15))
V2 <- VlnPlot(reference_list$CNON_SCZ, features = "percent.mt") + scale_y_continuous(limits = c(0, 15))
V1 + V2

# Run Neighbors and Clustering for SEPscL310
SEPscL310_seuratobj <- FindNeighbors(SEPscL310_seuratobj, reduction = "pca", dims = 1:30)
SEPscL310_seuratobj <- FindClusters(SEPscL310_seuratobj, resolution = 0.5)

# DimPlot
p1 <- DimPlot(SEPscL310_seuratobj, group.by = "seurat_clusters", label = TRUE) + NoLegend() + xlim(-15, 15) + ylim(-15, 20)

# Find Markers
SEPscL310_seuratobj.markers <- FindAllMarkers(SEPscL310_seuratobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(SEPscL310_seuratobj.markers, file = "~/310cellcyclefilter.csv")

# Rename Clusters
new.cluster.ids <- c(rep("Others", 21), "Mesenchymal", "Others", "Others")
names(new.cluster.ids) <- levels(SEPscL310_seuratobj)
SEPscL310_seuratobj <- RenameIdents(SEPscL310_seuratobj, new.cluster.ids)
SEPscL310_seuratobj$CellType <- Idents(SEPscL310_seuratobj)

p2 <- DimPlot(SEPscL310_seuratobj, group.by = "seurat_clusters", label = TRUE) + NoLegend() + xlim(-15, 15) + ylim(-15, 20)
p1 + p2

# Feature Plots
VlnPlot(SEPscL310_seuratobj, features = c("S100A1", "PIP"))
FeaturePlot(SEPscL310_seuratobj, features = "IGHA1")

# Find Transfer Anchors
anchors <- lapply(reference_list, function(ref) {
  FindTransferAnchors(reference = SEPscL310_seuratobj, query = ref, query.assay = "RNA", normalization.method = "SCT", k.filter = NA, k.anchor = 10, reference.reduction = "pca", dims = 1:30)
})

# Map Query Data
reference_list <- lapply(seq_along(reference_list), function(i) {
  MapQuery(anchorset = anchors[[i]], query = reference_list[[i]], reference = SEPscL310_seuratobj, refdata = list(celltype.1 = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")
})

# Predictions and Violin Plots
predictions <- lapply(anchors, function(anchor) {
  TransferData(anchorset = anchor, refdata = SEPscL310_seuratobj$seurat_clusters, dims = 1:30)
})

for (i in seq_along(reference_list)) {
  reference_list[[i]] <- AddMetaData(reference_list[[i]], metadata = predictions[[i]])
}

VlnPlot(reference_list[[1]], features = "prediction.score.Neural.Progenitors", group.by = "orig.ident")
VlnPlot(reference_list[[2]], features = "prediction.score.Neural.Progenitors", group.by = "orig.ident")

# Define Markers
np.markers.036 <- c("TAGLN", "IGFBP7", "COL1A2", "CALD1", "TPM2", "COL1A1", "COL3A1", "TPM1", "LGALS1", "SERPINB3", "VWF", "GAPDH", "ACTB", "FUT4", "CD24", "ITGB1")
np.markers.181 <- c("TAGLN", "IGFBP7", "COL1A2", "CALD1", "TPM2", "COL1A1", "COL3A1", "TPM1", "LGALS1", "VWF", "SNTN", "GAPDH", "ACTB", "DMBT1", "FUT4", "CD24", "ITGB1", "KRT5")
facet <- c("TAGLN", "IGFBP7", "COL1A2", "CALD1", "TPM2", "COL1A1", "COL3A1", "TPM1", "LGALS1")
hk.markers <- c("GAPDH", "ACTB")
np.markers <- c("TAGLN", "IGFBP7", "COL1A2", "CALD1", "TPM2", "COL1A1", "COL3A1", "TPM1", "LGALS1")
Features <- c("MC" = "TAGLN", "MC" = "IGFBP7", "MC" = "COL1A2", "MC" = "COL1A1", "MC" = "CALD1", "MC" = "TPM2", "MC" = "COL3A1", "MC" = "TPM1", "MC" = "LGALS1", "Housekeeping" = "GAPDH", "Housekeeping" = "ACTB")
messen.markers <- c("mesenchymal" = "CD44", "mesenchymal" = "VIM", "mesenchymal" = "NT5E", "mesenchymal" = "THY1", "mesenchymal" = "LUM", "mesenchymal" = "ANPEP", "mesenchymal" = "ALCAM", "mesenchymal" = "HSPA8")
messen.markers <- c("other mesenchymal markers" = "CD44", "other mesenchymal markers" = "VIM", "Must be expressed" = "NT5E", "Must be expressed" = "THY1", "Must be expressed" = "ENG", "other mesenchymal markers" = "ANPEP", "other mesenchymal markers" = "ALCAM", "other mesenchymal markers"
                    