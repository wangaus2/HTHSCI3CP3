#install.packages("here")

# Set environment

library("here")
library("dplyr")
setwd(here())

path = "./data"
raw_path = "./data/raw"
raw_file_paths = list.files(raw_path)
meta_path = "./data/metadata"
annotatedIDs = paste0(meta_path, "/", list.files(meta_path))

dat = sapply(
  raw_file_paths, 
  function(x) 
    read.delim(paste0(raw_path, '/', x), header = F))

dat = do.call('cbind', dat)

cat("[OK] Setting environment")

rownames(dat) <- dat[,1]
dat <- dat[, dat[1, ] != "Geneid"]
names(dat) <- substr(names(dat), 1, 3)

#removing 1st Row/Last Col
dat<- dat[-1,]
dat<- dat[,-18]

# converting to numeric
dat[] = lapply(dat, as.numeric)

#Handling Missing Values 
if (any(is.na(dat))) {
  cat("[FAIL] Missing entries detected")
} else {
  cat("[OK] No missing entries")
}

if (all(dat %% 1 == 0)) {
  cat("[OK] All data in Integers")
} else {
  cat("[FAIL] Not all data in Integers")
}

cat("[OK] Cleaning Data")

#Make annotations
geneID = read.delim(annotatedIDs, header = T)
annot = data.frame(sample_name = colnames(dat))
annot = cbind(annot, condition = c(rep('Control', 8), rep('Treatment', 9)))
rownames(annot) <- annot[, 1]
cat("[OK] Making Annotations")

#Write Files
out_dir <- file.path(path, "processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# write processed data
write.table(
  dat,
  file.path(out_dir, "dat_processed.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)

# write sample annotation
write.table(
  annot,
  file.path(out_dir, "sample_annotation.txt"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
)
  
cat(
  "[OK] Complete, outputs written to:",
  out_dir,
  "\n"
)