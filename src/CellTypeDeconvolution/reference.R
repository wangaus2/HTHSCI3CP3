source(here("src/CellTypeDeconvolution/reference_header.R"))

paths <- list.files(here("data/raw/PDAC_scRNA"), full.names = TRUE)
expression_matrices <- lapply(paths, Read10X)

file_names <- list.files(here("data/raw/PDAC_scRNA"))
seurat_objects <- CreateSeuratObjects(expression_matrices = expression_matrices, paths = file_names)


seurat_objects_filtered <- lapply(seurat_objects, qc)
