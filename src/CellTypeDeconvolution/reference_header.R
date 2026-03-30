qc <- function(seurat_object){
  seurat_object@active.assay = 'RNA'
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") 
  seurat_object <- subset(
    seurat_object, 
    subset = nFeature_RNA > 200 & 
      nFeature_RNA < 2500 & 
      percent.mt < 5
  )
  return(seurat_object)
}

CreateSeuratObjects <- function(expression_matrices, paths){
  seurat_objects <- list()
  for (i in seq_along(expression_matrices)){
    seurat_objects[[i]] <- CreateSeuratObject(counts = expression_matrices[i], 
                                              project = paths[i])
  }
  return(seurat_objects)
}
