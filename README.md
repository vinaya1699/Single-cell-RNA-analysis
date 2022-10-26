# Single-cell-RNA-analysis
# Libraries used in an analysis :- 
library(Seurat)
library(patchwork)
library(dplyr)

# Input Data
pbmc.data = Read10X(data.dir = "C:/Users/91973/Documents/filtered_gene_bc_matrices/hg19/")

pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc

pbmc.data[1:50, 1:10]

# Calculating PercentageFeatureSet
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

# Vlnplot of pbmc data
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
![vlnplot](https://user-images.githubusercontent.com/110582335/197968954-dace2fa5-f1d6-4f3f-afe4-0fe3b4b72a9a.png)

# Creating a subset of pbmc whereis nFeature_RNA is between 200 to 2500
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
pbmc = NormalizeData(pbmc)

# Find Variable Features of pbmc data
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data 
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

pbmc@assays$RNA@scale.data[1:50, 1:5]

# Run PCA using variable features of pbmc data
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find K Neighghbor of pbmc data
pbmc = FindNeighbors(pbmc, dims = 1:10)

# Plotting a dimheatmap
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Find Clusters in PBMC data
pbmc = FindClusters(pbmc, resolution = 0.5)

head(pbmc@meta.data)

# Use method UMAP 
pbmc = RunUMAP(pbmc, dims = 1:10)

# Plot a Dimplot of calculated UMAP
DimPlot(pbmc, reduction = "umap", label = T)
![dimplot2](https://user-images.githubusercontent.com/110582335/197971858-a7f48fbe-7c1f-46a6-a89b-e3add25daa6e.png)

# To find all markers
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(pbmc.markers)

# Assign Markers to pbmc clustered data

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a

# Print genes 
genes = a %>% pull(gene)
genes

#Plot a feature plot	
FeaturePlot(pbmc, features = genes[1:2])
![featureplot](https://user-images.githubusercontent.com/110582335/197972458-bee10a3b-9cb4-4ede-bc28-9b9a7af43fd1.png)

