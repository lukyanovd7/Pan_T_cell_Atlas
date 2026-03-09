library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(future)
library(qs)

source("/home/osuchalko/DeepSingleCell/sc_funcs.R")
source("/home/osuchalko/DeepSingleCell/harmony_runs/Integrate_layers_big_dataset-HARMONY.R")


print("------- Reading Data ---------")
samples <- qread("/projects/retroelements/abu/runs/rds_files_deepsc/cleaned_harmony_INTEGRATED_d30_res03_no7_no10.rds")
print("------- Postprocessing ---------")
samples <- seurat_postprocessing(samples, d=30, red="harmony", res=0.3)
qsave(samples, "/projects/retroelements/abu/runs/rds_files_deepsc/cleaned_harmony_POSTPROCESSED_d30_res03_no7_no10_2.rds")
