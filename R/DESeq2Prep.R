library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
  BiocManager::install("DESeq2", force = TRUE)
  
  BiocManager::install("apeglm")
  library(apeglm)

BiocManager::install("vsn")
library(vsn)
library(DESeq2)
library(here)

path <-here()

processed_data_path <- paste0(path, "/data/processed")
annot <- read.delim(paste0(processed_data_path, "/sample_annotation.txt")) %>% as.data.frame()
annot <- annot[, -1, drop = FALSE]
dat <- read.delim(paste0(processed_data_path, "/dat_processed.tsv"))


deseq2_obj <- DESeqDataSetFromMatrix(countData = dat,
                                     colData = annot,
                                     design = ~ condition)

minimum_sample_size = 3
keep = rowSums(counts(deseq2_obj) >= 10) >= minimum_sample_size
deseq2_obj = deseq2_obj[keep, ]

dds = DESeq(deseq2_obj)
res = results(dds, contrast = c('condition', 'Control', 'Treatment'))

summary(res, alpha = 0.05)
head(res)

vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds) #simple log transformation by log2(count + 1)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) #variance stabilized counts

library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))


library(dplyr)
library(ggplot2)
# Bar plot
bar_dat = data.frame(Direction = c('Upregulated', 'Downregulated'),
                     NSig = c(subset(res, padj < 0.05 & log2FoldChange > 0) %>% nrow(), 
                              subset(res, padj < 0.05 & log2FoldChange < 0) %>% nrow()))

ggplot(bar_dat, aes(x = Direction, y = NSig, fill = Direction)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(y = 'Number of Differentially Expressed Genes (FDR < 0.05)', x = 'Direction')

# Volcano Plot
res2 = as.data.frame(res) %>% 
  mutate(gene_symbol = as.integer(rownames(.)), 
         log10P = -log10(pvalue), 
         Direction = ifelse(log2FoldChange < 0 , 'Negative', 'Positive'), 
         FDR = ifelse(padj < 0.05, '< 0.05', '> 0.05'), 
         highlight = ifelse(padj < sort(padj)[31], 1, 0))
res2 = left_join(res2, geneID[, c('ENTREZID', 'SYMBOL')], by = c('gene_symbol' = 'ENTREZID'))
ggplot(subset(res2, !is.na(log10P) & !is.na(padj)), aes(x = log2FoldChange, y = log10P, fill = Direction, shape = FDR, 
                                                        alpha = FDR)) +
  geom_point() +
  scale_shape_manual(values = c('< 0.05' = 21,
                                '> 0.05' = 21)) +
  scale_alpha_manual(values = c('< 0.05' = 1, 
                                '> 0.05' = 0.25)) +
  scale_fill_manual(values = c('Negative' = "steelblue3",
                               'Positive' = "coral2")) +
  geom_hline(yintercept = subset(res2, padj > 0.05)$log10P %>% max(), linetype = 'dashed') +
  labs(x = expression(paste("Log"[2], "Fold Change")), y = expression(paste("-log"[10], "P")), 
       fill = expression('Direction')) +
  facet_grid(cols = vars('Responder vs. Non-Responder Differential Gene Expression')) +
  ggrepel::geom_label_repel( data=subset(res2, highlight == 1), aes(label=SYMBOL), size = 3,
                             max.overlaps = 120, direction = 'both', alpha = 0.8, color = 'black') +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14, face = 'bold'),
        strip.background = element_rect(fill = 'seashell2'),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 14), 
        legend.position = 'bottom')
