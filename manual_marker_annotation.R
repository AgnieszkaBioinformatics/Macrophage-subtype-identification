rm(list = ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(jackstraw)
library(sva)
library(Matrix)



## reading in the data

base_dir <- "C:/Users/aurin/Downloads/GSE136001_RAW"
samples <- c("GSM4039243_f-tumor-1", "GSM4039244_f-tumor-2", "GSM4039241_f-ctrl-1", "GSM4039242_f-ctrl-2", "GSM4039247_m-tumor-1", "GSM4039248_m-tumor-2", "GSM4039245_m-ctrl-1", "GSM4039246_m-ctrl-2")

for (sample in samples) {
  mtx_dir <- file.path(base_dir, paste0(sample, "-filtered-matrix.mtx"), 
                       paste0(sample, "-filtered-matrix.mtx"))
  features_dir <- file.path(base_dir, paste0(sample, "-filtered-features.tsv"), 
                            paste0(sample, "-filtered-features.tsv"))
  bar_dir <- file.path(base_dir, paste0(sample, "-filtered-barcodes.tsv"))
  
  data_mtx <- ReadMtx(mtx = mtx_dir, features = features_dir, cells = bar_dir)
  
  name <- strsplit(sample, "_")[[1]][2]
  name <- gsub("-", "_", name)
  assign(name, CreateSeuratObject(counts = data_mtx))
}



## merging Seurat objects (it's not integration step yet) to check QC simultaneously
merged <- merge(f_ctrl_1, y = c(f_ctrl_2, f_tumor_1, f_tumor_2, m_ctrl_1, m_ctrl_2, m_tumor_1, m_tumor_2),
                add.cell.ids = c("female_c1", "female_c2", "female_t1", "female_t2", "male_c1",
                                 "male_c2", "male_t1", "male_t2"))
merged$sample <- rownames(merged@meta.data)
merged@meta.data <- separate(merged@meta.data, col = "sample", into = c("sex", "type", "barcode"),
                             sep = "_")

f_ctrl <- merge(f_ctrl_1, y=f_ctrl_2, add.cell.ids = c("female_c1", "female_c2"))
m_ctrl <- merge(m_ctrl_1, y = m_ctrl_2, add.cell.ids = c("male_c1", "male_c2"))
f_tumor <- merge(f_tumor_1, y = f_tumor_2, add.cell.ids = c("female_t1", "female_t2"))
m_tumor <- merge(m_tumor_1, y = m_tumor_2, add.cell.ids = c("male_t1", "male_t2"))


## QC
merged$percent_mt <- PercentageFeatureSet(merged, pattern = "^mt-")
merged$percent_mt
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot1
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

filtered <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt < 5)


f_ctrl$percent_mt <- PercentageFeatureSet(f_ctrl, pattern = "^mt-")
filtered_f_ctrl <- subset(f_ctrl, 
                          subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt <5)

m_ctrl$percent_mt <- PercentageFeatureSet(m_ctrl, pattern = "^mt-")
filtered_m_ctrl <- subset(m_ctrl, 
                          subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt <5)

f_tumor$percent_mt <- PercentageFeatureSet(f_tumor, pattern = "^mt-")
filtered_f_tumor <- subset(f_ctrl, 
                          subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt <5)

m_tumor$percent_mt <- PercentageFeatureSet(m_tumor, pattern = "^mt-")
filtered_m_tumor <- subset(m_tumor, 
                           subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt <5)


## Normalization
normalized_log <- NormalizeData(filtered, normalization.method = "LogNormalize", scale.factor = 10000)

f_ctrl_n <- NormalizeData(filtered_f_ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
m_ctrl_n <- NormalizeData(filtered_m_ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
f_tumor_n <- NormalizeData(filtered_f_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
m_tumor_n <- NormalizeData(filtered_m_tumor, normalization.method = "LogNormalize", scale.factor = 10000)


## Feature selection
glioma_log <- FindVariableFeatures(normalized_log, selection.method = "vst", nfeatures = 2000)
top10_log <- head(VariableFeatures(glioma_log), 10)
plot_log1 <- VariableFeaturePlot(glioma_log)
plot_log2 <- LabelPoints(plot = plot_log1, points = top10_log, repel = TRUE)
plot_log2

f_ctrl_n <- FindVariableFeatures(f_ctrl_n, selection.method = "vst", nfeatures = 2000)
m_ctrl_n <- FindVariableFeatures(m_ctrl_n, selection.method = "vst", nfeatures = 2000)
f_tumor_n <- FindVariableFeatures(f_tumor_n, selection.method = "vst", nfeatures = 2000)
m_tumor_n <- FindVariableFeatures(m_tumor_n, selection.method = "vst", nfeatures = 2000)


## Scaling
all.genes_log <- rownames(glioma_log)
glioma_log <- ScaleData(glioma_log, features = all.genes_log)

all.f_ctrl_n <- rownames(f_ctrl_n)
f_ctrl_n <- ScaleData(f_ctrl_n, features = all.f_ctrl_n)
all.m_ctrl_n <- rownames(m_ctrl_n)
m_ctrl_n <- ScaleData(m_ctrl_n, features = all.m_ctrl_n)

all.f_tumor_n <- rownames(f_tumor_n)
f_tumor_n <- ScaleData(f_tumor_n, features = all.f_tumor_n)
all.m_tumor_n <- rownames(m_tumor_n)
m_tumor_n <- ScaleData(m_tumor_n, features = all.m_tumor_n)


##### before batch correction ####

## PCA
glioma_log <- RunPCA(glioma_log, features = VariableFeatures(object = glioma_log))
ElbowPlot(glioma_log)
DimPlot(glioma_log, reduction = "pca", split.by = "type") 
DimPlot(glioma_log, reduction = "pca", split.by = "sex") 
DimPlot(glioma_log, reduction = "pca", group.by = "sex")
DimPlot(glioma_log, reduction = "pca", group.by = "type")

f_ctrl <- RunPCA(f_ctrl_n, features = VariableFeatures(object = f_ctrl_n))
ElbowPlot(f_ctrl)
m_ctrl <- RunPCA(m_ctrl_n, features = VariableFeatures(object = m_ctrl_n))
ElbowPlot(m_ctrl)
f_tumor <- RunPCA(f_tumor_n, features = VariableFeatures(object = f_tumor_n))
ElbowPlot(f_tumor)
m_tumor <- RunPCA(m_tumor_n, features = VariableFeatures(object = m_tumor_n))
ElbowPlot(m_tumor)


js <- JackStraw(glioma_log, num.replicate = 100)
js <- ScoreJackStraw(js, dims = 1:20)
JackStrawPlot(js, dims = 1:20)


## KNN clustering
knn <- FindNeighbors(glioma_log, dims = 1:20)
knn <- FindClusters(knn, resolution = 0.5)


knn_f_ctrl <- FindNeighbors(f_ctrl, dims = 1:20)
knn_f_ctrl <- FindClusters(knn_f_ctrl, resolution = 0.5)

knn_m_ctrl <- FindNeighbors(m_ctrl, dims = 1:20)
knn_m_ctrl <- FindClusters(knn_m_ctrl, resolution = 0.5)

knn_f_tumor <- FindNeighbors(f_tumor, dims = 1:20)
knn_f_tumor <- FindClusters(knn_f_tumor, resolution = 0.5)

knn_m_tumor <- FindNeighbors(m_tumor, dims = 1:20)
knn_m_tumor <- FindClusters(knn_m_tumor, resolution = 0.8)


## kmeans clustering
embeddings <- Embeddings(m_ctrl, reduction = "pca")
km_mc <- kmeans(embeddings, centers = 12)

## UMAP 

umap <- RunUMAP(knn, dims = 1:20)
p1 <-DimPlot(umap, reduction = "umap", group.by = "sex")
p2 <- DimPlot(umap, reduction = "umap", group.by = "type")
p1+p2

DimPlot(umap, reduction = "umap")
        
umap_f_c <- RunUMAP(knn_f_ctrl, dims = 1:20)
umap_m_c <- RunUMAP(knn_m_ctrl, dims = 1:20)
umap_f_t <- RunUMAP(knn_f_tumor, dims = 1:20)
umap_m_t <- RunUMAP(knn_m_tumor, dims = 1:20)



## TSNE
tsne_knn <- RunTSNE(knn, dims = 1:20)
DimPlot(object=tsne_knn, reduction = 'tsne')

tsne_f_ctrl <- RunTSNE(knn_f_ctrl, dims = 1:20)
DimPlot(object=tsne_f_ctrl, reduction = 'tsne')
tsne_m_ctrl <- RunTSNE(knn_m_ctrl, dims = 1:20)
DimPlot(object=tsne_m_ctrl, reduction = 'tsne')

tsne_f_tumor <- RunTSNE(knn_f_tumor, dims = 1:20)
DimPlot(object=tsne_f_tumor, reduction = 'tsne')
tsne_m_tumor <- RunTSNE(knn_m_tumor, dims = 1:20)
DimPlot(object=tsne_m_tumor, reduction = 'tsne')


## BATCH CORRECTION
glioma_log <- IntegrateLayers(object = glioma_log, method = CCAIntegration, orig.reduction = "pca",
                              new.reduction = "integrated.cca", verbose = FALSE)

glioma_log[["RNA"]] <- JoinLayers(glioma_log[["RNA"]])



knn <- FindNeighbors(glioma_log, reduction = "integrated.cca", dims = 1:20)
knn <- FindClusters(knn, resolution = 1)


umap <- RunUMAP(knn, dims = 1:20, reduction = "integrated.cca")
DimPlot(umap, reduction = "umap", group.by = "type")
DimPlot(umap, reduction = "umap", group.by = "sex")
DimPlot(umap_m_t, reduction = "umap")


fc <- IntegrateLayers(object = umap_f_c, method = CCAIntegration, orig.reduction = "pca",
                      new.reduction = "integrated.cca", verbose = FALSE)

mc <- IntegrateLayers(object = umap_m_c, method = CCAIntegration, orig.reduction = "pca",
                      new.reduction = "integrated.cca", verbose = FALSE)

ft <- IntegrateLayers(object = umap_f_t, method = CCAIntegration, orig.reduction = "pca",
                      new.reduction = "integrated.cca", verbose = FALSE)

mt <- IntegrateLayers(object = umap_m_t, method = CCAIntegration, orig.reduction = "pca",
                      new.reduction = "integrated.cca", verbose = FALSE)

fc[["RNA"]] <- JoinLayers(fc[["RNA"]])
mc[["RNA"]] <- JoinLayers(mc[["RNA"]])
ft[["RNA"]] <- JoinLayers(ft[["RNA"]])
mt[["RNA"]] <- JoinLayers(mt[["RNA"]])

## Finding differentially expressed features (cluster biomarkers)
cluster2.markers <- FindAllMarkers(umap, only.pos = TRUE)

cluster.biomarkers <- data.frame(cluster2.markers)
write.table(cluster.biomarkers, file = "C:/Users/aurin/Downloads/GSE136001_results/cluster_biomarkers.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

VlnPlot(glioma_log, features = "Itga4", layer = "counts", log = TRUE, group.by = "sex")

cluster2.markers_fc <- FindAllMarkers(fc, only.pos = TRUE)
cluster2.markers_mc <- FindAllMarkers(mc, only.pos = TRUE)
cluster2.markers_ft <- FindAllMarkers(ft, only.pos = TRUE)
cluster2.markers_mt <- FindAllMarkers(mt, only.pos = TRUE)

## Marker based annotation

# homeostatic microglia
FeaturePlot(umap, features = c("Crybb1", "Cst3", "P2ry12", "Pros1"))

# mg2
FeaturePlot(umap, features = c("Jun", "Junb", "Jund", "Fos", "Egr1", "Klf6", 
                               "Aft3", "Nfkbia"))

# mg3
FeaturePlot(umap, features = c("Bmp2k", "Bhlhe41", "Ncoa3", "Notch2"))

# mg5
FeaturePlot(umap, features = c("Plp1", "Pltp", "Mbp"))

# mg6
FeaturePlot(umap, features = c("Cd63", "Cd9"))

# microglia
FeaturePlot(umap, features = c("Tmem119", "P2ry12", "Crybb1"))

# premature microglia
FeaturePlot(umap, features = c("Csf1", "Mcm5", "Ifit3", "Cst7", "Mif", "Ccl12", 
                               "Ccl3", "Ccl4"))

# disease associated microglia
FeaturePlot(umap, features = c("B2m", "H2-D1", "H2-K1", "H2-Oa", "H2-DMa", 
                               "Bst2", "Lgals3bp", "Ccl12" ))

# mg 8
FeaturePlot(umap, features = c("Stmn1", "Tubb5", "Tuba1b", "Cdk1", "Top2a"))

# monocytes
FeaturePlot(umap, features = c("Ly6c2", "Ccr2"))

 # NK
FeaturePlot(umap, features = c("Ncam1"))

# dendritic cells
FeaturePlot(umap, features = c("Dc", "Cd24a"))


# Mo/MΦ
FeaturePlot(umap, features = c("Ly6c2", "Ccr2", "Tgfbi"))


# macrophages
FeaturePlot(umap, features = c("Ly6c2", "Ifitm2", "Ifitm3", "S100a6"))

# t cells
FeaturePlot(umap, features = c("Cd3d", "Cd3e", "Cd3g"))

# cytotoxic t cell
FeaturePlot(umap, features = c("Cd8a", "Cd8b1"))

# b cell
FeaturePlot(umap, features = c("Cd19", "Ms4a1", "Sdc1"))

#bam Apoe, Ms4a7, and Mrc1
FeaturePlot(umap, features = c("Apoe", "Ms4a7", "Mrc1"))


fc_new_ids <- c("Mg4", "homeostatic microglia", "Mg2", "Mg3", "BAM", 
                "premature microglia", "6", "Mg6")

names(fc_new_ids) <- levels(umap_f_c)
umap_f_c_test <- RenameIdents(umap_f_c, fc_new_ids)
DimPlot(umap_f_c_test, reduction = "umap")

mc_new_ids <- c("homeostatic microglia", "homeostatic microglia", "Mg2", "Mg3", "BAM", 
                "monocytes", "NK", "dendritic cells", "premicroglia", "9")

names(mc_new_ids) <- levels(umap_m_c)
umap_m_c <- RenameIdents(umap_m_c, mc_new_ids)
DimPlot(umap_m_c, reduction = "umap")


ft_new_ids <- c("premature microglia", "homeostatic microglia", "Mg2", "Mg3", 
                "BAM", "disease associated microglia", "macrophages, Mo/MΦ", "NK")
names(ft_new_ids) <- levels(umap_f_t)
umap_f_t <- RenameIdents(umap_f_t, ft_new_ids)
DimPlot(umap_f_t, reduction = "umap")


mt_new_ids <- c("homeostatic microglia", "homeostatic microglia",
                "homeostatic microglia", "Mg2","Mg3", "Mg8", "BAM", 
                "premature microglia", "NK", "dendritic cells",
                "disease associated microglia", "macrophages, Mo/MΦ", 
                "monocytes", "13")
names(mt_new_ids) <- levels(umap_m_t)
umap_m_t <- RenameIdents(umap_m_t, mt_new_ids)
DimPlot(umap_m_t, reduction = "umap")


merged_new_ids <- c("Mg6", "Mg2", "homeostatic microglia", 
                    "disease associated microglia, premature microglia", "Mg3", "Mo/MΦ", "BAM",
                    "macrophages, monocytes", "T cells", "Mg8", "dendritic cells",
                    "11", "12", "13", "B cells", "cytotoxic T cells", "16")

names(merged_new_ids) <- levels(umap)
umap1 <- RenameIdents(umap, merged_new_ids)
DimPlot(umap1, reduction = "umap")

## preparing data for anndata convertion 


f_ctrl_ad <- GetAssayData(f_ctrl_2, assay = 'RNA', slot = 'counts')
writeMM(f_ctrl_ad, file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\fc_matrix.mtx")
write.table(data.frame('gene' = rownames(f_ctrl_ad)), 
            file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\fc_genenames.csv",
            quote = F, row.names = F, col.names = F
            )

f_ctrl_2$barcode <- colnames(f_ctrl_2)
write.csv(f_ctrl_2@meta.data, 
          file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\fc_metadata.csv",
          quote = F, row.names = F)


# m_ctrl
m_ctrl_ad <- GetAssayData(m_ctrl_1, assay = 'RNA', slot = 'counts')
writeMM(m_ctrl_ad, file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mc_matrix.mtx")
write.table(data.frame('gene' = rownames(m_ctrl_ad)), 
            file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mc_genenames.csv",
            quote = F, row.names = F, col.names = F
)

m_ctrl_1$barcode <- colnames(m_ctrl_1)
write.csv(m_ctrl_1@meta.data, 
          file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mc_metadata.csv",
          quote = F, row.names = F)

# ft
f_tumor_ad <- GetAssayData(f_tumor_1, assay = 'RNA', slot = 'counts')
writeMM(f_tumor_ad, file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\ft_matrix.mtx")
write.table(data.frame('gene' = rownames(f_tumor_ad)), 
            file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\ft_genenames.csv",
            quote = F, row.names = F, col.names = F
)

f_tumor_1$barcode <- colnames(f_tumor_1)
write.csv(f_tumor_1@meta.data, 
          file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\ft_metadata.csv",
          quote = F, row.names = F)

# mt
m_tumor_ad <- GetAssayData(m_tumor_1, assay = 'RNA', slot = 'counts')
writeMM(m_tumor_ad, file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mt_matrix.mtx")
write.table(data.frame('gene' = rownames(m_tumor_ad)), 
            file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mt_genenames.csv",
            quote = F, row.names = F, col.names = F
)

m_tumor_1$barcode <- colnames(m_tumor_1)
write.csv(m_tumor_1@meta.data, 
          file = "C:\\Users\\aurin\\Downloads\\GSE136001_results\\dane\\mt_metadata.csv",
          quote = F, row.names = F)


