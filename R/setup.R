if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", force=TRUE)
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install('clusterProfiler')
BiocManager::install("org.hs.eg.db")
cran_packages <- c(
  "dplyr",
  "pheatmap",
  "RColorBrewer",
  "here",
  "ggplot2"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

rm(cran_packages)
cat("setup complete")