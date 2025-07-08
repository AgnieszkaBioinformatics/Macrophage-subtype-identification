library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(celldex)
library(SingleR)
library(cellmarkeraccordion)
BiocManager::install("cellmarkeraccordion")

setwd("C:\\Users\\aurin\\Desktop\\Trento\\courses\\single cells")

load("QCC.rds")

####### QC
qcc <- CreateSeuratObject(counts = counts, 
                          project = "qcc", # name of the project
                          min.cells = 3,   # filter for genes (rows)
                          min.features = 200 # filter for cells (columns)
)

## acess the count matrix
layer_counts <- LayerData(qcc, assay = "RNA", layer = "counts")

## calculate mitochondrial genes
qcc[["percent_mt"]] <- PercentageFeatureSet(qcc, pattern = "^MT-")

VlnPlot(qcc, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, pt.size = 0.01)


plot1 <- FeatureScatter(qcc, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(qcc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

qcc <- 
  subset(qcc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 
         & percent_mt < 10) # 5%: 34 cells, 10%: 116 cells


## normalization
qcc <- 
  NormalizeData(qcc, normalization.method = "LogNormalize", scale.factor = 10000)


## feature selection
qcc <- FindVariableFeatures(qcc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(qcc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(qcc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot2



## cell cycle analysis
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
qcc <- CellCycleScoring(qcc, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
RidgePlot(qcc, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

## scaling
all.genes <- rownames(qcc)
qcc <- ScaleData(qcc, features = all.genes)


## PCA
qcc <- RunPCA(qcc, features = VariableFeatures(object = qcc))
ElbowPlot(qcc) # 20 PCs
DimHeatmap(qcc, dims = 20, cells = 30, balanced = TRUE)
qcc <- RunPCA(qcc, features = c(s_genes, g2m_genes))

# Visualize top genes associated with reduction components
VizDimLoadings(qcc, dims = 1:2, reduction = "pca")


DimPlot(qcc, reduction = "pca")



## clustering -> UMAP and tSNE not suited for small datasets
qcc <- FindNeighbors(qcc, dims = 1:20)
qcc <- FindClusters(qcc, resolution = 0.5)

qcc <- RunUMAP(qcc, dims = 1:20)
DimPlot(qcc, reduction = "umap")

set.seed(1)
qcc <- RunTSNE(qcc)
DimPlot(qcc, reduction = "tsne")

tsne_out <- data.frame(qcc@reductions$tsne@cell.embeddings)
#tsne_data <- data.table(tsne_out, layer_counts)



## kmeans
pca_emb <- Embeddings(qcc, "pca")[, 1:20] 
set.seed(123)
km <- kmeans(pca_emb, centers = 5) # optimal number
qcc$kmeans_clustering <- factor(km$cluster)
DimPlot(qcc, reduction = "tsne", group.by = "kmeans_clustering")

Idents(qcc) <- "kmeans_clustering"
## differentially expressed features


markers <- FindAllMarkers(qcc, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

markers %>%
  group_by(cluster) %>%
  slice_max(n=4,order_by=avg_log2FC)


## hierarchical clustering
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(qcc, features = top10$gene) + NoLegend()


## plot marker expression
VlnPlot(qcc, features = c("CCL5", "SERPINF1"),pt.size=0)
VlnPlot(qcc, features = c("CD2", "LYZ"),pt.size=0)
RidgePlot(qcc, features = c("LYZ"))
FeaturePlot(qcc, features = c("CCL5", "SERPINF1", "CD2", "LYZ"), order =T )
DotPlot(qcc, features = c("SPP1", "EGR2", "HBEGF",
                          "RGR", "RLBP1", "RPE65", "ITGB8",
                          "PXK", "CGAS", "AC007216.3", "SLC25A53",
                          "AC004540.1", "GUCA1A", "AP5S1", "ARL4D",
                          "CD3D", "CCL5", "IFITM1", "RPS27"))


#### annotation
ref_human <- celldex::HumanPrimaryCellAtlasData()
ref_mouse <- celldex::MouseRNAseqData()

View(data.frame(colData(ref_mouse)))
m <- data.frame(colData(ref_mouse))
length(unique(m$label.main)) # 18

h <- data.frame(colData(ref_human))
length(unique(h$label.main)) # 36

ann_human <- SingleR(test = layer_counts,
        ref = ref_human,
        labels = ref_human$label.main)

ann_mouse <- SingleR(test = layer_counts,
                     ref = ref_mouse,
                     labels = ref_mouse$label.main)

qcc$ann_human <- 
  ann_human$pruned.labels[match(rownames(qcc@meta.data), rownames(ann_human))]
DimPlot(qcc, reduction = "tsne", group.by = "ann_human")

qcc$ann_mouse <- 
  ann_mouse$pruned.labels[match(rownames(qcc@meta.data), rownames(ann_mouse))]
DimPlot(qcc, reduction = "tsne", group.by = "ann_mouse")

## assessing the quality of annotations
plotScoreHeatmap(ann_human)
plotScoreHeatmap(ann_mouse)


## using fine assignment
ann_human_f <- SingleR(test = layer_counts,
                     ref = ref_human,
                     labels = ref_human$label.fine)

ann_mouse_f <- SingleR(test = layer_counts,
                     ref = ref_mouse,
                     labels = ref_mouse$label.fine)

qcc$ann_human_f <- 
  ann_human_f$pruned.labels[match(rownames(qcc@meta.data), rownames(ann_human_f))]
DimPlot(qcc, reduction = "tsne", group.by = "ann_human_f")

qcc$ann_mouse_f <- 
  ann_mouse_f$pruned.labels[match(rownames(qcc@meta.data), rownames(ann_mouse_f))]
DimPlot(qcc, reduction = "tsne", group.by = "ann_mouse_f")


## assessing the quality of annotations
plotScoreHeatmap(ann_human_f)
plotScoreHeatmap(ann_mouse_f)

plotDeltaDistribution(ann_human_f)
plotDeltaDistribution(ann_human)

plotDeltaDistribution(ann_mouse_f)
plotDeltaDistribution(ann_mouse)



####### cellmarkeraccordion
