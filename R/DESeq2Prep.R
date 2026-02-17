library(here)
source(here("R", "header.R"))

#Compile DESeq2 Object
processed_data_path <- here("data/processed")
setwd(here("data"))
raw_file_paths <- list.files("raw")
annotatedIDs <- list.files("metadata/")
geneID <- read.delim(paste0("metadata/", annotatedIDs), header = T)

#ind <- grep("*M" , colnames(dat))
#dat <- dat[, -ind]
ids <- c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,8)

annot <- read.delim(paste0(processed_data_path, "/sample_annotation.txt")) %>% as.data.frame()
annot <- annot[, -1, drop = FALSE]
annot[,1] <- factor(annot[,1]) #Ensure Annotations are as factor
annot<- mutate(annot, cell_ids = ids)
annot[,2] <- factor(annot[,2]) #Ensure Annotations are as factor

#annot <- annot[-ind, ]
#rownames(annot)
#rownames(annot) <- c(1:11)

dat <- read.delim(paste0(processed_data_path, "/dat_processed.tsv"))


deseq2_obj <- DESeqDataSetFromMatrix(countData = dat,
                                     colData = annot,
                                     design = ~ cell_ids + condition)

minimum_sample_size = 2
keep = rowSums(counts(deseq2_obj) >= 20) >= minimum_sample_size
deseq2_obj = deseq2_obj[keep, ]


dds = DESeq(deseq2_obj)

res = results(dds, contrast = c('condition', 'Treatment', 'Control'))

summary(res, alpha = 0.05)
head(res)

dfres <- as.data.frame (res)
significant_values<- dfres[dfres$padj < 0.05 & !is.na(dfres$padj) & (dfres$log2FoldChange > 0.5 | dfres$log2FoldChange < -0.5) , ]

vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds) #simple log transformation by log2(count + 1)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) #variance stabilized counts
file <- here("data/processed/DEA_allsamples.csv")
write.csv(res, file, row.names = FALSE)

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
