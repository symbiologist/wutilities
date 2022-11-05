# Process and prepare 10x pbmc 10k dataset to use for testing functions
# Download raw data from 10x and extract to data directory
# https://www.10xgenomics.com/resources/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0

library(tidyverse)
library(Seurat)
library(here)

# Import
seuratobj <- CreateSeuratObject(counts = Read10X_h5(here('data/pbmc_10k_v3_filtered_feature_bc_matrix.h5')), 
                                min.cells = 100,
                                min.features = 1000)

# Add fake metadata for testing
seuratobj$Treatment <- c('Treated', 'Untreated') %>% sample(ncol(seuratobj), replace = TRUE)

# Standard pipeline
seuratobj <- seuratobj %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures()

# Dimensionality reduction
seuratobj <- seuratobj %>% RunPCA(npcs = 10) %>% RunUMAP(dims = 1:10)

# Clustering
seuratobj <- seuratobj %>% FindNeighbors() %>% FindClusters(algorithm = '4', resolution = 0.2)

seuratobj %>% DimPlot()

# Slim down for saving
seuratobj <- seuratobj %>% 
  DietSeurat(counts = TRUE, 
             data = TRUE,
             scale.data = FALSE,
             features = NULL,
             assays = NULL,
             dimreducs = c('pca', 'umap'),
             graphs = c('RNA_nn', 'RNA_snn'),
             misc = TRUE)

# Export
seuratobj %>% write_rds(here('data/seuratobj.rds'), compress = 'bz')
