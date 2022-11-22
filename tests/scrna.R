library(Seurat)
library(tidyverse)
library(wutilities)
library(here)

seuratobj <- read_rds(here('data/seuratobj.rds'))

seurat_feature(seuratobj, label_size = 5, color_package = 'carto', color_palette = 'Prism')
