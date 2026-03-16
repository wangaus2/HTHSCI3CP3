if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioconductor_packages <- c(
  "AnnotationHub",
  "clusterProfiler",
  "AnnotationDbi",
  "GOSemSim",
  "vsn",
  "apeglm",
  "DESeq2",
  "org.Hs.eg.db"
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
for (pkg in bioconductor_packages) {
  if (!require(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

rm(list = c("cran_packages", "bioconductor_packages"))
message("Setup complete")
