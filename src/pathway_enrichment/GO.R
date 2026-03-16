#Setup
    library(here)
    source(here("src/setup/header.R"))
    setwd(here())
    resfilter <- read.csv("results/DESeq/tables/DESeqResFiltered.csv")
    res <- read.csv("results/DESeq/tables/DESeqRes.csv")
    gene_ids <- res$gene_id
     up <- subset(resfilter, logFC < 0)$gene_id
    down <- subset(resfilter, logFC > 0)$gene_id
    
    sig_list <- list(up = up, down = down)
    
    go_res = lapply(sig_list, function(i) {
      
      res = lapply(c('BP', 'MF' ,'CC'), function(x) {
        message(paste0('Analyzing ', x))
        enrichGO(i,
                 OrgDb = org.Hs.eg.db, #If human, use org.Hs.eg.db
                 keyType = 'ENTREZID',
                 ont = x,
                 universe = gene_ids,
                 readable = F) 
      }) %>% setNames(c('BP', 'MF', 'CC'))
      
    })
    go_sig = lapply(c('up', 'down'), function(y) {
      
      i = go_res[[y]]
      
      sig = lapply(c('BP', 'MF', 'CC'), function(x) {
        
        df = subset(i[[x]]@result, p.adjust <= 0.05) 
        #Save only if there are significant results
        if (nrow(df) > 1) write.table(df, paste0( "../results/GO/tables/", y, '_GO', x, '_sig.txt'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t') 
        return(df)
        
      }) %>% setNames(c('BP', 'MF', 'CC'))
      
    }) %>% setNames(c('up', 'down'))
    
    go_summary = data.frame(matrix(0, ncol = 3, nrow = 2))
    rownames(go_summary) = c('up', 'down')
    colnames(go_summary) = c('BP', 'MF', 'CC')
    for (y in c('up', 'down')) {
      i = go_sig[[y]]
      for (n in c('BP', 'MF', 'CC')) {
        x = i[[n]]
        message(paste0(nrow(x),' GO BP terms are over-represented among ', y, '-regulated genes'))
        go_summary[y, n] = nrow(x)
      }
    }
    go_summary
    write.table(go_summary, '../results/GO/tables/go_summary.txt', col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t') #Save results
    