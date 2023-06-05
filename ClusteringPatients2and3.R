# Using Satija Guided Clustering Tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow) 
# on Espinoza et al. 2023 scRNAseq data! 

library(dplyr)
library(Seurat)
library(patchwork)

# Creating Seurat objects

patient2.data <- Read10X(data.dir = "/Users/danaalkalali/Downloads/Amit2")
patient3.data <- Read10X(data.dir = "/Users/danaalkalali/Downloads/Amit3")

combined <- CreateSeuratObject(counts = patient2.data, project = "Amit2", min.cells = 3, min.features = 200)
combined

patient3 <- CreateSeuratObject(counts = patient3.data, project = "Amit3", min.cells = 3, min.features = 200)
patient3

# Merging patients 2 and 3 into a combined dataset! :D 
combined <- merge(patient2, y = c(patient3), add.cell.ids = c("patient2", "patient3"), project = "23combined")
combined
head(colnames(combined))
table(combined$orig.ident)

-----
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# QC Seurat 
# Need to change parameters to same as Espinoza paper 

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Data Normalization 

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

combined <- FindVariableFeatures(combined)
top10 <- head(VariableFeatures(combined), 10)
# print(top10)
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

# Data Scaling! 

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

# Linear Dimension Reduction 

combined <- RunPCA(combined, features = VariableFeatures(object = combined))

print(combined[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(combined, dims = 1:2, reduction = "pca")

DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE)

# Determining Dimensionality 

combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)
JackStrawPlot(combined, dims = 1:15)

#PCs are strong because you are using preprocessed data dum dum! 

# Clustering time!!! Woohoo

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

head(Idents(combined), 5)

# Now time for non-linear DR! *_* 

combined <- RunUMAP(combined, dims = 1:10)
DimPlot(combined, reduction = "umap")

#Example, top 5 genes of cluster 2: 
cluster2.markers <- FindMarkers(combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

#Example, how cluster 5 is differentiated from clusters 0 and 3:
cluster5.markers <- FindMarkers(combined, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(combined, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
print(cluster0.markers)


# I chose the top 2 genes from the results to graph in a violin plot
VlnPlot(combined, features = c("MALAT1", "CXCR4"))

combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(combined, features = top10$gene) + NoLegend()
