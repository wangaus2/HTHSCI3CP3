vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds) #simple log transformation by log2(count + 1)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) #variance stabilized counts

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Add sample names (these are the "actual names")
pcaData$sample <- colnames(vsd)


# Plot with labels
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw()

# Bar plot
bar_dat = data.frame(Direction = c('Upregulated', 'Downregulated'),
                     NSig = c(subset(res, padj < 0.05 & log2FoldChange > 0.5) %>% nrow(), 
                              subset(res, padj < 0.05 & log2FoldChange < -0.5) %>% nrow()))
bar_dat <- data.frame(
  Direction = c("Upregulated", "Downregulated"),
  NSig = c(
    subset(resfilter, logFC < 0) %>% nrow(),
    subset(resfilter, logFC > 0) %>% nrow()
  )
)

ggplot(bar_dat, aes(x = Direction, y = NSig, fill = Direction)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(y = 'Number of Differentially Expressed Genes (FDR < 0.05)', x = 'Direction') +
  scale_fill_manual(values = c("Downregulated" = "steelblue3",
                               "Upregulated" = "coral2"))

# Volcano Plot
res2 = as.data.frame(res) %>% 
  mutate(gene_symbol = as.integer(rownames(.)), 
         log10P = -log10(pvalue), 
         Direction = ifelse(log2FoldChange < 0 , 'Upregulated', 'Downregulated'), 
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
  scale_fill_manual(values = c('Upregulated' = "coral2",
                               'Downregulated' = "steelblue3")) +
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
