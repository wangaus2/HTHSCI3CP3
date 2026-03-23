bioconductor_packages <- c(
  "AnnotationHub",
  "clusterProfiler",
  "AnnotationDbi",
  "GOSemSim",
  "vsn",
  "apeglm",
  "DESeq2",
  "org.Hs.eg.db",
  "GSVA",
  "UCSCXenaTools"
)
cran_packages <- c(
  "dplyr",
  "pheatmap",
  "RColorBrewer",
  "here",
  "ggplot2",
  "msigdbr",
  "ggrepel"
)

invisible(lapply(c(bioconductor_packages, cran_packages), library, character.only = TRUE))
rm(list = c("bioconductor_packages", "cran_packages"))
setwd(here())
