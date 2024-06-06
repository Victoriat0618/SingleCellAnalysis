# Load libraries
library(Seurat)
library(data.table)
library(tidyverse)
library(sctransform)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(presto)
library(monocle3)

# Assuming data is already loaded and processed in Seurat object 'seurat_obj'

# Extract meta data from Seurat object
cell_annotations <- as.data.frame(seurat_obj@meta.data)

# Extract gene names from PCA feature loadings from Seurat object
gene_annotation <- data.frame(gene_short_name = rownames(seurat_obj@reductions[["pca"]]@feature.loadings))
rownames(gene_annotation) <- rownames(seurat_obj@reductions[["pca"]]@feature.loadings)

# Assuming expression_matrix_t is correctly transposed with genes as rows and samples as columns
gene_annotation_corrected <- data.frame(gene_short_name = rownames(expression_matrix_t))
rownames(gene_annotation_corrected) <- rownames(expression_matrix_t)

# Load the info as cell data set into Monocle
cds <- new_cell_data_set(expression_matrix_t, cell_metadata = cell_annotations, gene_metadata = gene_annotation_corrected)

# Recreate partition
recreate.partition <- rep(1, length(cds@colData@rownames))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

# Assign cluster information
list_cluster <- seurat_obj@meta.data$seurat_clusters
names(list_cluster) <- colnames(seurat_obj)
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

# Extract UMAP coordinates from Seurat object and add to Monocle cell_data_set
umap_coordinates <- Embeddings(seurat_obj, "umap")
reducedDims(cds)[["UMAP"]] <- umap_coordinates

# Extract PCA feature loadings from Seurat
pca_feature_loadings <- seurat_obj@reductions[["pca"]]@feature.loadings

# Plot top genes contributing to PC1
loadings_df <- as.data.frame(pca_feature_loadings[, 1])
loadings_df$gene <- rownames(loadings_df)
colnames(loadings_df)[1] <- "loading"
ggplot(loadings_df, aes(x=reorder(gene, loading), y=loading)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Top Genes Contributing to PC1", x="Gene", y="Loading") +
  coord_flip()

# Run trajectory analysis
cds <- learn_graph(cds, use_partition = TRUE)

# Visualize data
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, group_label_size = 5) + 
  theme(legend.position = "right")
plot_cells(cds, color_cells_by = "cluster", label_branch_points = T, label_roots = T, group_label_size = 5)

# Order cells by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 5]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T, label_branch_points = T)

# Extract and visualize pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

# Differential expression in trajectory analysis
deg <- graph_test(cds, neighbor_graph = "principal_graph", cores = 32)
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
deg_ids <- row.names(subset(deg, q_value < 0.05))

# Group trajectory genes into modules and visualize heatmap
gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")

# Define genes of interest
WntGenes <- c("FZD2", "FZD7", "ROR2", "LPR5", "LPR6", "SFRP1","SFRP2", "CTNNB1","Wnt5A", "Wnt5B", "Notch2",
              "JAG1", "PSEN1" ,"PSEN2", "ADAMS17", "PSENEN", "APH1A")

# Visualize gene expression
RidgePlot(seurat_obj, features = WntGenes, sort = TRUE) + scale_x_continuous(c(0, 5))
FeaturePlot(seurat_obj, features = WntGenes)

# Visualize trajectory of certain genes
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("FZD2", "FZD7", "ROR2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime")

# Plot pseudotime
plot_cells(cds, genes = WntGenes, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)

# Manually adjust certain branches
cds_sub <- choose_graph_segments(cds)
