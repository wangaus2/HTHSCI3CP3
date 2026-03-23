library(here)
source(here("src/setup/header.R"))

data("XenaData")

cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% 
  XenaScan("TCGA Pancreatic Cancer")

cli_query = cohort %>% # download clinical data
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>% 
  XenaDownload(destdir = here("data/raw"))# download to current directory

ge = cohort %>% # download gene expression data
  filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq") %>%
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>% 
  XenaDownload(destdir = here("data/raw"))


cli = XenaPrepare(
  subset(
    cli_query, grepl('PAAD', datasets)
  )
)

exp.data <- XenaPrepare(subset(ge, grepl('PAAD', datasets)))
exp.matrix <- as.matrix(exp.data[,-1])
rownames(exp.matrix) <- exp.data$sample

#rm (list = c("ge", "XenaData", "cli", "cli_query"))

## Load in Test Set
res = read.csv(here("results/DESeq/tables/DESeqResFiltered.csv"))

res$gene_id <- as.character(res$gene_id)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = res$gene_id,
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multiVals = "first")

## test sets will be up and dn set from DE seq analysis

res.up <- res[res$logFC < 0 , "symbol"]
res.up <- res.up[!is.na(res.up)]
res.dn <- res[res$logFC > 0, "symbol"]
res.dn <- res.dn[!is.na(res.dn)]

gs <- list(
  down = res.dn,
  up = res.up
)

## Performing GSVA analysis
gsva_obj = gsvaParam(
  exprData = exp.matrix, 
  geneSets = gs, 
  #minSize = 5, maxSize = 500, 
  kcdf = 'Gaussian'
  )
gsva_h = gsva(gsva_obj) %>% t() %>% as.data.frame()
gsva_h$sampleID = rownames(gsva_h)


## Clinical data
clinical = cli$PAAD_clinicalMatrix
primary = subset(clinical, sample_type == 'Primary Tumor')

#Matching sampleIDs between clinical data matrix and gsva result matrix
gsva_h <- gsva_h[match(primary$sampleID, gsva_h$sampleID),]
gsva_h <- gsva_h[-which(is.na(gsva_h$sampleID)),] ## removing mismatched values

primary = primary[match(gsva_h$sampleID, primary$sampleID), ] ## removing mismatched values


gsva_h = cbind(gsva_h,
               hpg = primary$neoplasm_histologic_grade)
gsva_h$hpg = factor(gsva_h$hpg, levels = c('G1', 'G2', 'G3', 'G4'))
    missingvalues <- which(is.na(gsva_h$hpg))
    gsva_h <- gsva_h[-missingvalues, ]
    primary <- primary[-missingvalues, ]

anova_up = aov(up ~ hpg, data = gsva_h)
summary(anova_up)
anova_dn = aov(down ~ hpg, data = gsva_h)
summary(anova_dn)


## Cox Proportional Hazard Regression
survival = cli$PAAD_survival.txt
survival = survival[match(gsva_h$sampleID, survival$sample), ]
s = Surv(time = survival$PFI.time, event = survival$PFI)
reg = coxph(s ~ up + down, data = gsva_h)
summary(reg)


## Multivariate Regression
gsva_h = cbind(gsva_h, 
               age = primary$age_at_initial_pathologic_diagnosis, 
               sex = primary$gender)

m_reg = coxph(s ~ up + age + sex, data = gsva_h)
summary(m_reg)
