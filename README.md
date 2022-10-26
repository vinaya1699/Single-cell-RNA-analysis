# Single-cell-RNA-analysis
# Libraries used in an analysis :- 
library(Seurat)
library(patchwork)
library(dplyr)

# Input Data
pbmc.data = Read10X(data.dir = "C:/Users/91973/Documents/filtered_gene_bc_matrices/hg19/")

pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc

#Output of pbmc 
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)

pbmc.data[1:50, 1:10]


#Output of pbmc.data[1:50, 1:10]
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

# Calculating PercentageFeatureSet
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

#Output of head(pbmc@meta.data)
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

# Creating a subset of pbmc whereis nFeature_RNA is between 200 to 2500
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
pbmc = NormalizeData(pbmc)

#Output : Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|

# Find Variable Features of pbmc data
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Outpput : 
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|

# Scale data 
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

pbmc@assays$RNA@scale.data[1:50, 1:5]

#Output :- 
AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 AAACCGTGTATGCG-1
AL627309.1         -0.05812316      -0.05812316      -0.05812316      -0.05812316      -0.05812316
AP006222.2         -0.03357571      -0.03357571      -0.03357571      -0.03357571      -0.03357571
RP11-206L10.2      -0.04166819      -0.04166819      -0.04166819      -0.04166819      -0.04166819
RP11-206L10.9      -0.03364562      -0.03364562      -0.03364562      -0.03364562      -0.03364562
LINC00115          -0.08223981      -0.08223981      -0.08223981      -0.08223981      -0.08223981
NOC2L              -0.31717081      -0.31717081      -0.31717081      -0.31717081      -0.31717081
KLHL17             -0.05344722      -0.05344722      -0.05344722      -0.05344722      -0.05344722
PLEKHN1            -0.05082183      -0.05082183      -0.05082183      -0.05082183      -0.05082183
RP11-54O7.17       -0.03308805      -0.03308805      -0.03308805      -0.03308805      -0.03308805
HES4               -0.23376818      -0.23376818      -0.23376818      -0.23376818      -0.23376818
RP11-54O7.11       -0.03768905      -0.03768905      -0.03768905      -0.03768905      -0.03768905
ISG15              -0.83530282      -0.83530282       0.39223510       2.21976210      -0.83530282

# Run PCA using variable features of pbmc data
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Output :- 
PC_ 1 
Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
	   FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
	   PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
	   CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
	   MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
PC_ 2 
Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
	   HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
	   BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
	   CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
	   TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC

# Find K Neighghbor of pbmc data
pbmc = FindNeighbors(pbmc, dims = 1:10)

# Plotting a dimheatmap
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Find Clusters in PBMC data
pbmc = FindClusters(pbmc, resolution = 0.5)

#Output :- 
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 2638
Number of edges: 95965

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8723
Number of communities: 9
Elapsed time: 0 seconds

head(pbmc@meta.data)
## head(pbmc@meta.data)
orig.ident nCount_RNA nFeature_RNA percent.mt RNA_snn_res.0.5 seurat_clusters
AAACATACAACCAC-1 SeuratProject       2419          779  3.0177759               2               2
AAACATTGAGCTAC-1 SeuratProject       4903         1352  3.7935958               3               3
AAACATTGATCAGC-1 SeuratProject       3147         1129  0.8897363               2               2
AAACCGTGCTTCCG-1 SeuratProject       2639          960  1.7430845               1               1
AAACCGTGTATGCG-1 SeuratProject        980          521  1.2244898               6               6
AAACGCACTGGTAC-1 SeuratProject       2163          781  1.6643551               2               2

# Use method UMAP 
pbmc = RunUMAP(pbmc, dims = 1:10)

#Output :- 
13:42:23 UMAP embedding parameters a = 0.9922 b = 1.112
13:42:23 Read 2638 rows and found 10 numeric columns
13:42:23 Using Annoy for neighbor search, n_neighbors = 30
13:42:23 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
13:42:23 Writing NN index file to temp file C:\Users\91973\AppData\Local\Temp\Rtmps57ttS\file2be41307d25
13:42:23 Searching Annoy index using 1 thread, search_k = 3000
13:42:24 Annoy recall = 100%
13:42:25 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
13:42:25 Initializing from normalized Laplacian + noise (using irlba)
13:42:25 Commencing optimization for 500 epochs, with 105124 positive edges
Using method 'umap'
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
13:42:32 Optimization finished

# Plot a Dimplot of calculated UMAP
DimPlot(pbmc, reduction = "umap", label = T)
![dimplot2](https://user-images.githubusercontent.com/110582335/197971858-a7f48fbe-7c1f-46a6-a89b-e3add25daa6e.png)

# To find all markers
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(pbmc.markers)

#head(pbmc.markers)
              p_val avg_log2FC pct.1 pct.2     p_val_adj cluster  gene
RPS12 1.806317e-144  0.7350248 1.000 0.991 2.477183e-140       0 RPS12
RPS6  7.135900e-142  0.6798622 1.000 0.995 9.786173e-138       0  RPS6
RPS27 5.257820e-140  0.7207819 0.999 0.992 7.210575e-136       0 RPS27
RPL32 4.229582e-136  0.6115515 0.999 0.995 5.800448e-132       0 RPL32
RPS14 1.799019e-130  0.6199183 1.000 0.994 2.467175e-126       0 RPS14
RPS25 5.507298e-123  0.7442491 0.997 0.975 7.552709e-119       0 RPS25

# Assign Markers to pbmc clustered data

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a
#a
#A tibble: 18 x 7
#Groups:   cluster [9]
       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
 1 1.74e-109       1.07 0.897 0.593 2.39e-105 0       LDHB    
 2 1.17e- 83       1.33 0.435 0.108 1.60e- 79 0       CCR7    
 3 0               5.57 0.996 0.215 0         1       S100A9  
 4 0               5.48 0.975 0.121 0         1       S100A8  
 5 7.99e- 87       1.28 0.981 0.644 1.10e- 82 2       LTB     
 6 2.61e- 59       1.24 0.424 0.111 3.58e- 55 2       AQP3    
 7 0               4.31 0.936 0.041 0         3       CD79A   
 8 9.48e-271       3.59 0.622 0.022 1.30e-266 3       TCL1A   
 9 1.17e-178       2.97 0.957 0.241 1.60e-174 4       CCL5    
10 4.93e-169       3.01 0.595 0.056 6.76e-165 4       GZMK    
11 3.51e-184       3.31 0.975 0.134 4.82e-180 5       FCGR3A  
12 2.03e-125       3.09 1     0.315 2.78e-121 5       LST1    
13 1.05e-265       4.89 0.986 0.071 1.44e-261 6       GZMB    
14 6.82e-175       4.92 0.958 0.135 9.36e-171 6       GNLY    
15 1.48e-220       3.87 0.812 0.011 2.03e-216 7       FCER1A  
16 1.67e- 21       2.87 1     0.513 2.28e- 17 7       HLA-DPB1
17 7.73e-200       7.24 1     0.01  1.06e-195 8       PF4     
18 3.68e-110       8.58 1     0.024 5.05e-106 8       PPBP  

# Print genes 
genes = a %>% pull(gene)
genes
#genes
 [1] "LDHB"     "CCR7"     "S100A9"   "S100A8"   "LTB"      "AQP3"     "CD79A"    "TCL1A"    "CCL5"     "GZMK"    
[11] "FCGR3A"   "LST1"     "GZMB"     "GNLY"     "FCER1A"   "HLA-DPB1" "PF4"      "PPBP"    

#Plot a feature plot	
FeaturePlot(pbmc, features = genes[1:2])
![featureplot](https://user-images.githubusercontent.com/110582335/197972458-bee10a3b-9cb4-4ede-bc28-9b9a7af43fd1.png)

