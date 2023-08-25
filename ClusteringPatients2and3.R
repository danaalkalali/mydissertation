# Using Satija Guided Clustering Tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow) 
# on Espinoza et al. 2023 scRNAseq data! 

require(knitr)
install.packages("sctransform")

library(dplyr)
library(Seurat)
library(patchwork)
library(knitr)

# Creating Seurat objects

patient2.data <- Read10X(data.dir = "/Users/danaalkalali/Downloads/Amit2")
patient3.data <- Read10X(data.dir = "/Users/danaalkalali/Downloads/Amit3")

patient2 <- CreateSeuratObject(counts = patient2.data, project = "Amit2", min.cells = 3, min.features = 200)
patient2

patient3 <- CreateSeuratObject(counts = patient3.data, project = "Amit3", min.cells = 3, min.features = 200)
patient3

# Merging patients 2 and 3 into a combined dataset! :D 

combined <- merge(patient2, y = c(patient3), add.cell.ids = c("patient2", "patient3"), project = "23combined")
combined
head(colnames(combined))
table(combined$orig.ident)


#####

# Pre-processing: QC and cell selection (or the new term I have coined: cell-ection)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

head(combined@meta.data, 5)

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Viewing feature-feature relationships :D 

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Data Normalization 

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(combined), 10)
top10

# print(top10)
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

#####

# Data Scaling before the big girl stuff 

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

# Linear Dimension Reduction!!!!!!!!!!!!!!!!!

combined <- RunPCA(combined, features = VariableFeatures(object = combined))

print(combined[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(combined, dims = 1:3, reduction = "pca")

DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1:2, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE)

# Determining Dimensionality 

combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)
JackStrawPlot(combined, dims = 1:20) #since can only do up to 20 


# Clustering time!!! Woohoo. Changed resolution to 1.3 to match Espinoza and PCs to 35) - still haven't tried it out 

combined <- FindNeighbors(combined, dims = 1:35)
combined <- FindClusters(combined, resolution = 1.3)

head(Idents(combined), 5)
# What do those levels in the results ^ 0-8 represent? 

# Now time for non-linear DR! *_* 

combined <- RunUMAP(combined, dims = 1:35) #should it be 35 like above?
DimPlot(combined, reduction = "umap")

# Top 10 genes of the clusters, however I'm using FindAllMarkers function instead of 
# FindMarkers compared to the tutorial :D: 

cluster0.markers <- FindMarkers(combined, ident.1 = 0, min.pct = 0.25)
print(paste(rownames(x = cluster0.markers)[1:50], collapse = ", "))

cluster1.markers <- FindMarkers(combined, ident.1 = 1, min.pct = 0.25)
print(paste(rownames(x = cluster1.markers)[1:50], collapse = ", "))

cluster2.markers <- FindMarkers(combined, ident.1 = 2, min.pct = 0.25)
print(paste(rownames(x = cluster2.markers)[1:50], collapse = ", "))

cluster3.markers <- FindMarkers(combined, ident.1 = 3, min.pct = 0.25)
print(paste(rownames(x = cluster3.markers)[1:50], collapse = ", "))

cluster4.markers <- FindMarkers(combined, ident.1 = 4, min.pct = 0.25)
print(paste(rownames(x = cluster4.markers)[1:50], collapse = ", "))

cluster5.markers <- FindMarkers(combined, ident.1 = 5, min.pct = 0.25)
print(paste(rownames(x = cluster5.markers)[1:50], collapse = ", "))

cluster6.markers <- FindMarkers(combined, ident.1 = 6, min.pct = 0.25)
print(paste(rownames(x = cluster6.markers)[1:50], collapse = ", "))

cluster7.markers <- FindMarkers(combined, ident.1 = 7, min.pct = 0.25)
print(paste(rownames(x = cluster7.markers)[1:50], collapse = ", "))

cluster8.markers <- FindMarkers(combined, ident.1 = 8, min.pct = 0.25)
print(paste(rownames(x = cluster8.markers)[1:50], collapse = ", "))

cluster9.markers <- FindMarkers(combined, ident.1 = 9, min.pct = 0.25)
print(paste(rownames(x = cluster9.markers)[1:50], collapse = ", "))

cluster10.markers <- FindMarkers(combined, ident.1 = 10, min.pct = 0.25)
print(paste(rownames(x = cluster10.markers)[1:50], collapse = ", "))

cluster11.markers <- FindMarkers(combined, ident.1 = 11, min.pct = 0.25)
print(paste(rownames(x = cluster11.markers)[1:50], collapse = ", "))

cluster12.markers <- FindMarkers(combined, ident.1 = 12, min.pct = 0.25)
print(paste(rownames(x = cluster12.markers)[1:50], collapse = ", "))

cluster13.markers <- FindMarkers(combined, ident.1 = 13, min.pct = 0.25)
print(paste(rownames(x = cluster13.markers)[1:50], collapse = ", "))

cluster14.markers <- FindMarkers(combined, ident.1 = 14, min.pct = 0.25)
print(paste(rownames(x = cluster14.markers)[1:50], collapse = ", "))

cluster15.markers <- FindMarkers(combined, ident.1 = 15, min.pct = 0.25)
print(paste(rownames(x = cluster15.markers)[1:50], collapse = ", "))

cluster16.markers <- FindMarkers(combined, ident.1 = 16, min.pct = 0.25)
print(paste(rownames(x = cluster16.markers)[1:50], collapse = ", "))

cluster17.markers <- FindMarkers(combined, ident.1 = 17, min.pct = 0.25)
print(paste(rownames(x = cluster17.markers)[1:50], collapse = ", "))

cluster18.markers <- FindMarkers(combined, ident.1 = 18, min.pct = 0.25)
print(paste(rownames(x = cluster18.markers)[1:50], collapse = ", "))


#Other fun stuff, 

#How cluster 5 is differentiated from clusters 0 and 3:
cluster5.markers <- FindMarkers(combined, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

#Markers for every cluster compared to all remaining cells, reporting only positive
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# ROC test shows power of classification of the marker from 0 to 1 where 0 = random and 1 = perfect 
cluster0.markers <- FindMarkers(combined, ident.1 = 1, logfc.threshold = 0.3, test.use = "roc", only.pos = TRUE)
print(cluster0.markers)

FeaturePlot(combined, features = c("CXCR5", "IL6", "FCRL4", "CD27", "IL10", "CD70"))
FeaturePlot(combined, features = c("LTB"))
            
VlnPlot(combined, features = c("LTB"))

combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(combined, features = top10$gene) + NoLegend()


