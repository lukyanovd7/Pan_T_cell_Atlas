# The pan-T-cell axis: overlapping transcriptional identities across CD4⁺, CD8⁺, γδ, and MAIT T cells

## The code can be used to reproduce integration, postprocessing and figures from the manuscript

Data Format
Single-cell objects are currently saved in the .qs format using the qs R package. 
Existing .rds files can be re-saved via qsave(). 
The .qs format offers substantially faster read/write performance compared to .rds, so it is recommended in this case.

sc_funcs.R
Contains functions for both processing and post-processing of Seurat objects. 
Processing is applied to individual objects prior to integration, whereas post-processing is performed on the already-integrated object.

Integrate_layers_big_dataset-HARMONY.R
Performs multi-sample integration using the Harmony algorithm.

r_process_slurm.R
Executes integration and postprocessing, assuming that all samples have been saved as Seurat objects within a single .qs file. 
Adjusts TCR-seq metadata format, spliting cdr3 for two chains.

r_remove...
Implements re-integration, i.e. the recommended for comparative analysis with added/removed data (samples, subsets, etc.). 
The workflow reads the integrated object, splits it into individual samples by the RNA Assay, and repeats the integration pipeline from scratch. 
This is the script you should use if you intend to insert a custom preprocessing or comparison step prior to processing. 
Note that some objects contain CITE-seq ADT Assay, while others do not. 

