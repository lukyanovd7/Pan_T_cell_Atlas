library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(future)
library(DoubletFinder)
library(qs)

source("/home/osuchalko/DeepSingleCell/sc_funcs.R")
source("/home/osuchalko/DeepSingleCell/Integrate_layers_big_dataset-HARMONY.R")
plan("multicore", workers = 30)
options(future.globals.maxSize = 500 * 1024^3)


print("------- Reading Data ---------")
samples <- qread("/projects/single_cell/deeeepSC/objects/cleaned_processed_d30_res_to2.qs")
samples@meta.data$old_clusters <- samples$RNA_snn_res.2


print("------- Subset and split ---------")
samples[["RNA"]] <- split(samples[["RNA"]], f = samples$sample_id)


print("------- Preprocessing ---------")
samples <- seurat_processing(samples, norm = FALSE)


print("------- Integrating Samples ---------")
samples.integrated <- integrate_harmony(samples)
qsave(samples.integrated, "/projects/single_cell/deeeepSC/objects/remove_cd4_cd8/integrated_harmony.qs")


# print("------- Read ---------")
# samples.integrated <- qread("/projects/single_cell/deeeepSC/objects/cleaned_integrated_harmony.qs")

# print("Remove BAD")    
# samples.integrated@meta.data[is.na(samples.integrated$orig.ident),]$orig.ident <- "BAD"
# samples.integrated <- subset(samples.integrated, subset = orig.ident != "BAD")

print("------- Postprocess ---------")
samples.integrated <- seurat_postprocessing(samples.integrated, d=30, red="harmony")

print("SAVE")
qsave(samples.integrated, "/projects/single_cell/deeeepSC/objects/remove_cd4_cd8/processed_d30_res_to2.qs")
