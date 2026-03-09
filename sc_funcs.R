resize <- function(w,h) {
    options(repr.plot.width = w, repr.plot.height = h)
}

seurat_processing <- function(so, norm = T){

    print("-------Preprocessing-------")
    if (norm == T){
        print("-------Normalizing-------")
        so <- NormalizeData(so)
    }
    so <- FindVariableFeatures(so)
    so <- ScaleData(so)
    var_genes = VariableFeatures(object = so)
    
    print("-------Excluding TR V,D,J and other genes from variable genes-------")
    unwanted_genes <- "^TR[ABDG.][VDJC.][0-9]+|^TRAB[VDJC.][0-9]+|^TRGD[VDJC.][0-9]+|^MT-|^HLA|^IG[HKL.]|XIST"
    # unwanted_genes <- "^TR[ABDG.][VDJC.][0-9]+|^TRAB[VDJC.][0-9]+|^TRGD[VDJC.][0-9]+|^MT-|^HLA|^IG[HKL.]|XIST|^CD4$|^CD8[AB.]$"
    var_genes = var_genes[-grep(unwanted_genes, var_genes, perl = TRUE)]
    VariableFeatures(so) <- var_genes
    
    print("-------Running PCA-------")
    so <- RunPCA(so, features = var_genes)
    return (so)
}


seurat_postprocessing <- function(so, d, red) {
    so <- FindNeighbors(
        so,
        reduction = red, 
        dims=1:d
    )
    so <- FindClusters(
        so, 
        resolution=seq(1.2, 2, by=0.1)
    )
    so <- RunUMAP(
        so, 
        dims = 1:d, 
        reduction = red
    )
    all.genes <- rownames(so)
    so <- ScaleData(so, features = all.genes)
    return (so)
}


immune_features <- function(so) {
    p <- FeaturePlot(so,
                     features = c("CD3D", #  Tcells
                                  "CD4",
                                  "CD8A",
                                  "IL7R",
                                  "FOXP3",
                                  "MKI67",
                                  "KLRF1", #KLRF1 - NK
                                  "GZMB", "GZMA", "GNLY"
                                  # "MS4A1", #Bcells
                                  # "CD14", #CD14 - monocyte
                                  # "FCGR3A", #CD16 - monocyte and NK,
                                  # "HBA1" #HBA1 - erytro
                                 ),
                     # cols = c("lightgrey", "red"),
                     raster = TRUE,
                     label = TRUE,
                     label.size = 10,
                     repel = TRUE,
                     combine = FALSE
                    )
    return (p)
}


immune_dot <- function(so) {
    p <- DotPlot(so,
                 features = c('CD3D', 'CD3E', 
                         "CD4",
                         'CD8A', 'CD8B', 
                         'GNLY', 'NKG7', "GZMB", "GZMA",
                         "CD1D", "ZBTB16",
                         'IL7R', 'CCR6', 'CCR7',
                         'IL23R', 'IL17', 'MKI67',
                         # 'CD79A', 'MS4A1', "CD19",
                         'FCGR3A', 'KLRB1', "NCAM1", # NK
                         # 'LYZ', 'CD14', # CD14 Mono
                         'CST3', "IL3RA", # DC
                         "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2",
                         "TNF", "TIGIT", "IL17A")
                )+ RotatedAxis()
       # ) + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)    
    return (p)
}