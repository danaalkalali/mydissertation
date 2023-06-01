#Using Satija Guided Clustering Tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow) 
#on Espinoza et al. 2023 scRNAseq data! 

library(dplyr)
library(Seurat)
library(patchwork)

#Creating Seurat object 

patient2.data <- Read10X(data.dir = "/Users/danaalkalali/Downloads/Amit2")

patient2 <- CreateSeuratObject(counts = patient2.data, project = "Amit2", min.cells = 3, min.features = 200)
patient2

patient2[["percent.mt"]] <- PercentageFeatureSet(patient2, pattern = "^MT-")

#QC Seurat 
#Question: how will we choose # of features and mtDNA?

VlnPlot(patient2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(patient2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(patient2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

patient2 <- subset(patient2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Data Normalization 

patient2 <- NormalizeData(patient2, normalization.method = "LogNormalize", scale.factor = 10000)

patient2 <- FindVariableFeatures(patient2)
top10 <- head(VariableFeatures(patient2), 10)
#print(top10)
plot1 <- VariableFeaturePlot(patient2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

#Data Scaling! 

all.genes <- rownames(patient2)
patient2 <- ScaleData(patient2, features = all.genes)

#Linear Dimension Reduction 

patient2 <- RunPCA(patient2, features = VariableFeatures(object = patient2))

print(patient2[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(patient2, dims = 1:2, reduction = "pca")

DimPlot(patient2, reduction = "pca")
DimHeatmap(patient2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(patient2, dims = 1:15, cells = 500, balanced = TRUE)

#Determining Dimensionality 

patient2 <- JackStraw(patient2, num.replicate = 100)
patient2 <- ScoreJackStraw(patient2, dims = 1:20)
JackStrawPlot(patient2, dims = 1:15)

#All these PCs seem quite strong right?! (1-15) :) i.e. solid curve above dashed line. Should we add even more (compared to the tutorial?)

#Clustering time!!! Woohoo

patient2 <- FindNeighbors(patient2, dims = 1:10)
patient2 <- FindClusters(patient2, resolution = 0.5)

head(Idents(patient2), 5)
# What do those levels in the results ^ 0-8 represent? 

# Now time for non-linear DR! *_* 

patient2 <- RunUMAP(patient2, dims = 1:10)
DimPlot(patient2, reduction = "umap")

#Example, top 5 genes of cluster 2: 
cluster2.markers <- FindMarkers(patient2, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

#Example, how cluster 5 is differentiated from clusters 0 and 3:
#Question: do we want to keep the minimum percent 25% when comparing clusters? 
#Question: how many types of cluster comparisons do we want to make?
cluster5.markers <- FindMarkers(patient2, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


patient2.markers <- FindAllMarkers(patient2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
patient2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(patient2, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
print(cluster0.markers)
#Question: How do I interpret these results exactly?

#I chose the top 2 genes from the results to graph in a violin plot
VlnPlot(patient2, features = c("MALAT1", "CXCR4"))

patient2.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(patient2, features = top10$gene) + NoLegend()
#Getting an error ^^^ 

#Question: Does the last part of the tutorial apply to this dataset (assigning cell type identity to clusters)?
