library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

# helper function: convert one GO result to dataframe and keep top 10
prep_go_df <- function(go_obj, direction, ontology) {
  df <- as.data.frame(go_obj)
  
  if (nrow(df) == 0) {
    return(data.frame())
  }
  
  df %>%
    mutate(
      direction = direction,
      ontology = ontology,
      Description = str_wrap(Description, width = 40)
    ) %>%
    arrange(p.adjust) %>%              # top terms by adjusted p-value
    slice_head(n = 10)
}

# make the 6 dataframes
up_BP   <- prep_go_df(go_res$up$BP,   "KRAS-UP",   "BP")
up_MF   <- prep_go_df(go_res$up$MF,   "KRAS-UP",   "MF")
up_CC   <- prep_go_df(go_res$up$CC,   "KRAS-UP",   "CC")
down_BP <- prep_go_df(go_res$down$BP, "KRAS-DN", "BP")
down_MF <- prep_go_df(go_res$down$MF, "KRAS-DN", "MF")
down_CC <- prep_go_df(go_res$down$CC, "KRAS-DN", "CC")
max_fdr <- max(-log10(go_df$p.adjust), na.rm = TRUE)

common_scale <- scale_fill_gradient(
  name = expression(-log[10](FDR)),
  limits = c(0, max_fdr),
  low = "grey80",
  high = "red"
)

plot_go_bar <- function(df) {
  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = "Empty")
    )
  }
  df$Description <- stringr::str_wrap(df$Description, width = 18)
  ggplot(
    df,
    aes(
      x = FoldEnrichment,
      y = reorder(Description, FoldEnrichment),
      fill = -log10(p.adjust)
    )
  ) +
    geom_col() +
    common_scale +
    labs(
      title = paste(unique(df$direction), unique(df$ontology)),
      x = "Fold Enrichment",
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14)
    )
}

p1 <- plot_go_bar(up_BP)
p2 <- plot_go_bar(up_MF)
p3 <- plot_go_bar(up_CC)
p4 <- plot_go_bar(down_BP)
p5 <- plot_go_bar(down_MF)
p6 <- plot_go_bar(down_CC)

#combined_plot <- ((p1 | p2 | p3) / (p4 | p5 | p6)) +
  #plot_layout(guides = "collect")
combined_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "a",
    theme = theme(
      plot.title = element_text(hjust = 0, size = 16, face = "bold")
    )
  )

combined_plot & theme(
  legend.position = "right", 
  legend.title = element_text(size = 20),
  legend.text  = element_text(size = 14)
)
                      legend.text  = element_text(size = 14))