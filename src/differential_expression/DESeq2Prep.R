#Setup
    library(here)
    source(here("src/setup/header.R"))
    setwd(here("data"))

#Compile DESeq2 Object
    annotatedIDs = list.files(here("data/raw/metadata/"), full.names = TRUE)
    geneID <- read.delim(annotatedIDs, header = T)
    annot  <- read.delim(here("data/processed/sample_annotations.txt"))
    dat    <- read.delim(here("data/processed/counts.tsv"))
    annot$condition  <- factor(annot$condition, levels = c("Control", "Treatment"))
    annot$cell_line  <- factor(annot$cell_line)
    deseq2_obj <- DESeqDataSetFromMatrix(countData = dat,
                                         colData = annot,
                                         design = ~ cell_line + condition)
    minimum_sample_size = 2
    keep = rowSums(counts(deseq2_obj) >= 20) >= minimum_sample_size
    deseq2_obj = deseq2_obj[keep, ]
    
    dds = DESeq(deseq2_obj)
    
    res = as.data.frame(results(dds, contrast = c('condition', 'Treatment', 'Control')))
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id")
    colnames(res_df)[3] <- "logFC"
    colnames(res_df)[7] <- "FDR"
    res_df <-res
    res_df <- filter(res_df, !is.na(FDR) & !is.na(pvalue) & !is.na(stat))
    
    resfilter <- res_df %>%
      filter(abs(logFC) > 0.5, FDR < 0.05)

    
    outfiles <- c(here("results/DESeq/tables/DESeqRes.csv"), here("results/DESeq/tables/DESeqResFiltered.csv"))
    write.csv(res_df, file = outfiles[1], row.names = FALSE)
    write.csv(resfilter, file = outfiles[2], row.names = FALSE)
