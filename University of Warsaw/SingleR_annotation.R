library(SingleR)
library(celldex)
library(pheatmap)
library(Seurat)

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



## preprocess the data
# List of objects to process
objects <- list(f_ctrl_1, f_ctrl_2, m_ctrl_1, m_ctrl_2, f_tumor_1, f_tumor_2, m_tumor_1, m_tumor_2)

# Iterate over each object to perform the desired steps
for (i in 1:length(objects)) {
  # Add percent_mt column
  objects[[i]]$percent_mt <- PercentageFeatureSet(objects[[i]], pattern = "^mt-")
  
  # Subset the object based on the specified conditions
  objects[[i]] <- subset(objects[[i]], 
                         subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mt < 5)
  
  # Normalize the data
  objects[[i]] <- NormalizeData(objects[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", nfeatures = 2000)
  
  # Scale the data
  all_features <- rownames(objects[[i]])
  objects[[i]] <- ScaleData(objects[[i]], features = all_features)
  
  # Run PCA
  objects[[i]] <- RunPCA(objects[[i]], features = VariableFeatures(object = objects[[i]]))
  objects[[i]] <- RunUMAP(objects[[i]], dims = 1:20)
}

# Assign back to individual variables if needed
f_ctrl_1 <- objects[[1]]
f_ctrl_2 <- objects[[2]]
m_ctrl_1 <- objects[[3]]
m_ctrl_2 <- objects[[4]]
f_tumor_1 <- objects[[5]]
f_tumor_2 <- objects[[6]]
m_tumor_1 <- objects[[7]]
m_tumor_2 <- objects[[8]]


# Load the reference data
ref <- celldex::MouseRNAseqData()

# List of UMAP objects
umap_objects <- list(umap_f_c1, umap_f_c2, umap_m_c1, umap_m_c2, umap_f_t1, umap_f_t2, umap_m_t1, umap_m_t2)

# Process each UMAP object
for (i in 1:length(umap_objects)) {
  # Get assay data
  counts <- GetAssayData(umap_objects[[i]], slot = 'counts')
  
  # Run SingleR
  pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)
  
  # Assign labels
  umap_objects[[i]]$singleR.labels <- pred$labels[match(rownames(umap_objects[[i]]@meta.data), rownames(pred))]
  
  # Plot UMAP
  DimPlot(umap_objects[[i]], reduction = 'umap', group.by = 'singleR.labels')
}

# Assign back to individual variables if needed
umap_f_c1 <- umap_objects[[1]]
umap_f_c2 <- umap_objects[[2]]
umap_m_c1 <- umap_objects[[3]]
umap_m_c2 <- umap_objects[[4]]
umap_f_t1 <- umap_objects[[5]]
umap_f_t2 <- umap_objects[[6]]
umap_m_t1 <- umap_objects[[7]]
umap_m_t2 <- umap_objects[[8]]


