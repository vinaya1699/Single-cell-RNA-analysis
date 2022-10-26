# Single-cell-RNA-analysis
# Libraries used in an analysis :- 
library(Seurat)
library(patchwork)
library(dplyr)

# Input Data
pbmc.data = Read10X(data.dir = "C:/Users/91973/Documents/filtered_gene_bc_matrices/hg19/")

pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc

# Output of pbmc 
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)

pbmc.data[1:50, 1:10]
# Output of pbmc.data[1:50, 1:10]
50 x 10 sparse Matrix of class "dgCMatrix"
   [[ suppressing 10 column names ‘AAACATACAACCAC-1’, ‘AAACATTGAGCTAC-1’, ‘AAACATTGATCAGC-1’ ... ]]
                                 
MIR1302-10    . . . . . . . . . .
FAM138A       . . . . . . . . . .
OR4F5         . . . . . . . . . .
RP11-34P13.7  . . . . . . . . . .
RP11-34P13.8  . . . . . . . . . .
AL627309.1    . . . . . . . . . .
RP11-34P13.14 . . . . . . . . . .
RP11-34P13.9  . . . . . . . . . .
AP006222.2    . . . . . . . . . .
RP4-669L17.10 . . . . . . . . . .
OR4F29        . . . . . . . . . .
RP4-669L17.2  . . . . . . . . . .
RP5-857K21.15 . . . . . . . . . .
RP5-857K21.1  . . . . . . . . . .

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)
# Output of head(pbmc@meta.data)
                    orig.ident nCount_RNA nFeature_RNA percent.mt
AAACATACAACCAC-1 SeuratProject       2419          779  3.0177759
AAACATTGAGCTAC-1 SeuratProject       4903         1352  3.7935958
AAACATTGATCAGC-1 SeuratProject       3147         1129  0.8897363
AAACCGTGCTTCCG-1 SeuratProject       2639          960  1.7430845
AAACCGTGTATGCG-1 SeuratProject        980          521  1.2244898
AAACGCACTGGTAC-1 SeuratProject       2163          781  1.6643551

# Vlnplot of pbmc data
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
![vlnplot](https://user-images.githubusercontent.com/110582335/197968954-dace2fa5-f1d6-4f3f-afe4-0fe3b4b72a9a.png)

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc

pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

pbmc@assays$RNA@scale.data[1:50, 1:5]

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)
head(pbmc@meta.data)

pbmc = RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = T)

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(pbmc.markers)

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a

genes = a %>% pull(gene)
genes

FeaturePlot(pbmc, features = genes[1:2])

