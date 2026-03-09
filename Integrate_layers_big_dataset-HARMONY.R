integrate_harmony <- function(so) {
    all.merged <- so
    all.integrated <- IntegrateLayers(
          object = all.merged, method = HarmonyIntegration,
          orig.reduction = "pca", new.reduction = "harmony",
          verbose = FALSE
    )
    print("Integrate Layers done")
    all.integrated[["RNA"]] <- JoinLayers(all.integrated[["RNA"]])
    print("Join Layers done")
    # saveRDS(all.integrated, "/projects/retroelements/abu/runs/rds_files_deepsc/DEEP_INTEGRATED_HARMONY_JOINT.rds")
    return (all.integrated)
}
