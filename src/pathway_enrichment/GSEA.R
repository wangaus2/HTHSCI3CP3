#Setup
library(here)
source(here("src/setup/header.R"))
setwd(here())

#Retrieve Hallmark gene set
hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H" 
)
gene_sets <- split(hallmark$entrez_gene, hallmark$gs_name)

#Perform GSEA on our processed results
res <- read.csv("results/DESeq/tables/DESeqRes.csv")
res <- res[order(-res$stat),]
ranked_list <- res$stat
names(ranked_list) <- res$gene_id

#Perform GSEA on study processed results
res_study <- read.csv("data/processed/NIHMS2007906-supplement-Supplementary_Data_S1.csv")
res_study <- res_study[order(-res_study$logFC),]
res_study <- res_study[!is.na(res_study$entrezgene_id), ]
res_study <- res_study[!duplicated(res_study$entrezgene_id), ]
ranked_list_study <- res_study$logFC
names(ranked_list_study) <- res_study$entrezgene_id


gsea <- fgsea(pathways = gene_sets, stats = ranked_list)
gsea_study <-fgsea(pathways = gene_sets, stats = ranked_list_study)
merged <- merge(
  gsea[, c("pathway","NES","padj")],
  gsea_study[, c("pathway","NES","padj")],
  by = "pathway",
  suffixes = c("_klomp","_study")
)

gsea <- gsea[order(-gsea$NES)]
gsea <- gsea %>%
  mutate(log10p = -log10(padj))

gsea_study <- gsea_study[order(-gsea_study$NES)]
gsea_study <- gsea_study %>%
  mutate(log10p = -log10(padj))

label_setsDN <- c(
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS")
label_setsUP <-c(
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)


gsea <- gsea %>%
  mutate(
    idx = row_number(),
    pathway_label = gsub("^HALLMARK_", "", pathway)
  )

dn_df <- gsea %>% filter(pathway %in% label_setsDN)
up_df <- gsea %>% filter(pathway %in% label_setsUP)

gsea_study <- gsea_study %>%
  mutate(
    idx = row_number(),
    pathway_label = gsub("^HALLMARK_", "", pathway)
  )
dn_df_study <- gsea_study %>% filter(pathway %in% label_setsDN)
up_df_study <- gsea_study %>% filter(pathway %in% label_setsUP)

ggplot(gsea, aes(x = idx, y = NES)) +
  geom_point(aes(size = -log10(padj)), alpha = 0.5, colour = "blue2") +
  
  scale_size_continuous(
    name = expression(atop("Significance", -log[10](padj))),
    range = c(3, 11),
    breaks = c(0, 5, 10, 15, 20),
    limits = c(0, 40)
  )+
  geom_text_repel(
    data = dn_df,
    aes(label = pathway_label),
    direction = "x",
    hjust = 1,
    xlim = c(min(gsea$idx) +1, min(gsea$idx) + 35),
    segment.size = 0.4,
    segment.alpha = 0.8,
    size = 3
  ) + 

geom_text_repel(
  data = up_df,
  aes(label = pathway_label),
  direction = "x",
  hjust = 0,
  xlim = c(max(gsea$idx) - 0.5, max(gsea$idx) -40),
  segment.size = 0.4,
  segment.alpha = 0.8,
  size = 3
) + 
  theme_classic() +
  labs(
    x = "Hallmark gene sets",
    y = "Normalized Enrichment Score (NES)",
    size = expression(Significance~(-log[10](pval)))
  )
  