library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
# library(future)
library(qs)

source("/home/osuchalko/DeepSingleCell/sc_funcs.R")
source("/home/osuchalko/DeepSingleCell/Integrate_layers_big_dataset-HARMONY.R")

# plan("multicore", workers = 30)
# options(future.globals.maxSize = 200 * 1024^3)
# options(repr.plot.width = 15, repr.plot.height = 15)


# print("------- Reading Data ---------")
all_samples.list <- qread("/projects/single_cell/deeeepSC/objects/all_samples_list.qs")


# print("------- Merging Samples ---------")
all.merged <- merge(all_samples.list[[1]], 
                    y = all_samples.list[2:length(all_samples.list)],
                    project = "deeeeepSC")


# print("------- Preprocessing ---------")
all.merged <- seurat_processing(all.merged, norm = TRUE)
qsave(all.merged, "/projects/single_cell/deeeepSC/objects/check_results/all_samples_merged_preprocessed.qs")


# print("------- Integrating Samples ---------")
# all.merged <- qread("/projects/single_cell/deeeepSC/objects/all_samples_merged_preprocessed.qs")
all.integrated <- integrate_harmony(all.merged)
qsave(all.integrated, "/projects/single_cell/deeeepSC/objects/check_results/all_samples_integrated_harmony.qs")


print("------- Postprocessing ---------")
# all.integrated <- qread("/projects/single_cell/deeeepSC/objects/all_samples_integrated_harmony.qs")
all.integrated <- seurat_postprocessing(all.integrated, d=30, red="harmony")
all.integrated@meta.data <- all.integrated@meta.data %>% separate(CTaa, c('TCR1aa', 'TCR2aa'), sep = "_", remove = FALSE)
all.integrated@meta.data <- all.integrated@meta.data %>% separate(CTnt, c('TCR1nt', 'TCR2nt'), sep = "_", remove = FALSE)
qsave(all.integrated, "/projects/single_cell/deeeepSC/objects/check_results/all_samples_processed_d30resolutions.qs")
