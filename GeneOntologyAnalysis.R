library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)

# Read data
df <- read.csv("~/ConditioncomparisonASPC.csv", header = TRUE)
original_gene_list <- df$avg_log2FC
names(original_gene_list) <- df$X

# Prepare gene list for GSEA
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# Gene set enrichment analysis (GSEA)
GO_results <- gseGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  pAdjustMethod = "none"
)

# Plot GSEA results
gseaplot(GO_results, by = "all", title = GO_results$Description[17], geneSetID = 17)

# Dotplot for GSEA results
dotplot(GO_results, showCategory = 20, split = ".sign") + facet_grid(.~.sign)

# Save GSEA results to CSV
write.csv(GO_results@result, file = "GOASPClist.csv")

####################### Over-representation test #######################

# Read data
sigs <- read.csv("~/ConditioncomparisonASPC.csv", header = TRUE)
row.names(sigs) <- sigs$X

# Filter genes with log fold change > 0.5
genes_to_test <- rownames(sigs[sigs$avg_log2FC > 0.5,])

# Over-representation test
enrichGO_results <- enrichGO(
  gene = genes_to_test,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

# Plot over-representation results
goplot(enrichGO_results)
barplot(enrichGO_results, drop = TRUE, showCategory = 20, title = "GO Biological Pathways", font.size = 8)
dotplot(enrichGO_results, showCategory = 20)

# Compare clusters (optional, if you have multiple gene lists)
GL <- list(sigs$ID)
compare <- compareCluster(GL, fun = "enrichGO", OrgDb = org.Mm.eg.db)

# Additional GSEA plot
gseaplot(enrichGO_results, by = "all", title = enrichGO_results$Description[790], geneSetID = 1)
Key Points:
  Loading Required Packages: The script checks and installs the required packages (BiocManager, clusterProfiler, AnnotationDbi, org.Mm.eg.db) if they are not already installed.
Reading and Preparing Data: Data is read from a CSV file, and a gene list is prepared for GSEA by removing NAs and sorting in decreasing order.
GSEA: The gseGO function is used to perform GSEA, and results are visualized using gseaplot and dotplot.
Saving Results: GSEA results are saved to a CSV file.
Over-representation Test: Genes with a log fold change greater than 0.5 are filtered, and the enrichGO function is used to perform the over-representation test. Results are visualized using goplot, barplot, and dotplot.
Compare Clusters: An optional section is included for comparing clusters using the compareCluster function.
Additional Plot: An additional GSEA plot is created for a specific gene set.
This code structure ensures that each step is clearly defined and easy to follow, making it more maintainable and understandable for future use.

make comment more clear
ChatGPT
Sure, here is the revised version with clearer comments for better understanding:
  
  r
Copy code
# Load required packages
# Check if BiocManager is installed and install it if not
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")

# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)

# Read in the data from a CSV file
df <- read.csv("~/ConditioncomparisonASPC.csv", header = TRUE)

# Extract the average log2 fold change values
original_gene_list <- df$avg_log2FC

# Set the gene names as the names of the values
names(original_gene_list) <- df$X

# Remove NA values from the gene list
gene_list <- na.omit(original_gene_list)

# Sort the gene list in decreasing order, required for clusterProfiler
gene_list <- sort(gene_list, decreasing = TRUE)

# Perform Gene Set Enrichment Analysis (GSEA) using GO biological process terms
GO_results <- gseGO(
  gene = gene_list,        # Gene list for GSEA
  OrgDb = org.Mm.eg.db,    # Organism database
  keyType = "SYMBOL",      # Gene identifier type
  ont = "BP",              # Ontology (Biological Process)
  minGSSize = 3,           # Minimum gene set size
  maxGSSize = 800,         # Maximum gene set size
  pvalueCutoff = 0.05,     # P-value cutoff
  verbose = TRUE,          # Print messages
  pAdjustMethod = "none"   # P-value adjustment method
)

# Plot the GSEA results for a specific gene set
gseaplot(GO_results, by = "all", title = GO_results$Description[17], geneSetID = 17)

# Plot the dotplot for the top 20 categories
dotplot(GO_results, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign)

# Save the GSEA results to a CSV file
write.csv(GO_results@result, file = "GOASPClist.csv")

####################### Over-representation test #######################

# Read in the data from a CSV file
sigs <- read.csv("~/ConditioncomparisonASPC.csv", header = TRUE)
row.names(sigs) <- sigs$X

# Filter for genes with log2 fold change > 0.5
genes_to_test <- rownames(sigs[sigs$avg_log2FC > 0.5, ])

# Perform over-representation test using GO biological process terms
enrichGO_results <- enrichGO(
  gene = genes_to_test,    # Genes to test
  OrgDb = org.Mm.eg.db,    # Organism database
  keyType = "SYMBOL",      # Gene identifier type
  ont = "BP"               # Ontology (Biological Process)
)

# Plot the over-representation results
goplot(enrichGO_results)

# Plot a barplot of the top 20 categories
barplot(enrichGO_results, drop = TRUE, showCategory = 20, title = "GO Biological Pathways", font.size = 8)

# Plot a dotplot of the over-representation results
dotplot(enrichGO_results, showCategory = 20)

# Compare clusters (if you have multiple gene lists)
GL <- list(sigs$ID)
compare <- compareCluster(GL, fun = "enrichGO", OrgDb = org.Mm.eg.db)

# Additional GSEA plot for a specific gene set
gseaplot(enrichGO_results, by = "all", title = enrichGO_results$Description[790], geneSetID = 1)