library(Seurat)
library(tidyverse)
library(wutilities)
library(here)
library(harmony)

seuratobj <- read_rds(here('data/seuratobj.rds'))

seed <- 123

test1 <- seuratobj %>% FindNeighbors() %>% FindClusters(resolution = 0.8) %>% RunUMAP(dims = 1:10, seed.use = seed)
test1 %>% seurat_feature()
test1$seurat_clusters %>% table()
test1@reductions$umap@cell.embeddings %>% head()

###
test <- seuratobj %>% ScaleData()

test <- test %>% RunHarmony(dims = 1:10, group.by.vars = 'Treatment', assay.use = 'RNA')

test@reductions$harmony@cell.embeddings %>% head()

test <- test %>% FindNeighbors(reduction = 'harmony', dims = 1:10)

test <- test %>% RunUMAP(reduction = 'harmony', dims = 1:10)

test %>% seurat_feature()