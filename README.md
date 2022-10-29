# Single-cell-RNA-analysis
# Libraries used in an analysis :- 
library(Seurat)
library(patchwork)
library(dplyr)

# Input Data
pbmc_data = Read10X(data.dir = "C:/Users/91973/Documents/filtered_gene/rms/")

pbmc = CreateSeuratObject(counts = pbmc_data, min.cells = 3, min.features = 200)
pbmc

pbmc.data[1:50, 1:10]

# forming set of all genes starting with MT- as a set of mitochondrial genes
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
#The number of unique genes and total molecules are automatically calculated and are stored in the object meta data.
head(pbmc@meta.data)   

# Vlnplot of pbmc data
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
![vlnplot](https://user-images.githubusercontent.com/110582335/198816889-a9f44260-d1e7-419e-bd9d-34f3845b0d48.png)
#nFeature_RNA is the number of genes detected in each cell. 
#nCount_RNA is the total number of gene detected within a cell. 
#And each dot in the following plots represents a cell

# To Find relationships between features such as percent.mt,nCount_RNA,nFeature_RNA in here.
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
![featureplot1](https://user-images.githubusercontent.com/110582335/198821164-0b5a3229-9288-40ef-8bf8-d41faed5c00b.png)

plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
![featureplot2](https://user-images.githubusercontent.com/110582335/198821171-e3db500d-6730-4446-aa09-bbbe0180d36d.png)

#The number 0 and 0.94 above each plot denotes the correlations between x-axis and y-axis

# Creating a subset of pbmc whereis nFeature_RNA is between 200 to 2500
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
pbmc = NormalizeData(pbmc)

# Find Variable Features of pbmc data
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# To identify the 10 most highly variable genes
top10 = head(VariableFeatures(pbmc), 10)
top10

# To plot variable features with and without labels
plot1 = VariableFeaturePlot(pbmc)
plot1
![plot_without_labels](https://user-images.githubusercontent.com/110582335/198821355-27b0a0dc-ba44-4518-b4c4-b32549da5bab.png)


plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
![plot_with_labels](https://user-images.githubusercontent.com/110582335/198821772-9eb55922-8296-46b0-85e6-3872d4ef8842.png)



# To scale data 
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

pbmc@assays$RNA@scale.data[1:50, 1:5]

# Run PCA using variable features of pbmc data
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find K Neighghbor of pbmc data
pbmc = FindNeighbors(pbmc, dims = 1:10)

# Plotting a dimheatmap
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
![dimPlot](https://user-images.githubusercontent.com/110582335/198816974-77e065ad-ad68-4b77-8aaa-0aa6943060e7.png)



# Find Clusters in PBMC data
pbmc = FindClusters(pbmc, resolution = 0.5)

head(pbmc@meta.data)

# Use method UMAP 
pbmc = RunUMAP(pbmc, dims = 1:10)

# Plot a Dimplot of calculated UMAP
DimPlot(pbmc, reduction = "umap", label = T)
![dimplot1](https://user-images.githubusercontent.com/110582335/198817020-50bc7ec2-aa61-41c7-bb27-3d76e5b12daf.png)

#Using UMAP method (K - N algorithm ) clusters are created which represent clusters of cells with features that are very similar to each other  

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
![dimplot2](https://user-images.githubusercontent.com/110582335/198817099-a36e892b-1967-426b-9153-780f2ff46530.png)
# Interpretation :-
MYO16 and KCTD16 are highly expressed genes in above single-cell RNA data .
