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

res.up2 <- res[order(res$logFC, decreasing = FALSE), "symbol"]
res.up2 <- res.up2[!is.na(res.up2)]
res.up2 <- res.up2[1:200]

res.dn <- res[res$logFC > 0, "symbol"]
res.dn <- res.dn[!is.na(res.dn)]

res.dn2 <- res[order(res$logFC, decreasing = TRUE), "symbol"]
res.dn2 <- res.dn2[!is.na(res.dn2)]
res.dn2 <- res.dn2[1:200]

gs <- list(
  down = res.dn2,
  up = res.up2
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

#histologic grade
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
library(survival)

survival = cli$PAAD_survival.txt
survival = survival[match(gsva_h$sampleID, survival$sample), ]
s = Surv(time = survival$PFI.time, event = survival$PFI)
reg = coxph(s ~ up + down, data = gsva_h)
summary(reg)



## Multivariate Regression
gsva_h = cbind(gsva_h, 
               age = primary$age_at_initial_pathologic_diagnosis, 
               sex = primary$gender,
               grade= primary$neoplasm_histologic_grade)

m_reg = coxph(s ~ up + age + sex + grade, data = gsva_h)
summary(m_reg)

#```{r Visualize}------------------------------------------------
library(ggbeeswarm)

#```{r Visualize ANOVA}

#Visualize difference in score between grades
library(ggpubr)
ggplot(subset(gsva_h, !is.na(hpg)), aes(x = hpg, y = up, color = hpg)) +
  geom_boxplot() +
  geom_quasirandom(alpha = 0.2) +
  labs(y = 'Gene Set Score', x = 'Neoplasm Histologic Grade', color = '') +
  theme_bw() +
  facet_grid(cols = vars('Difference in Gene Set Score Between Neoplasm Grades')) +
  theme(strip.text = element_text(face = 2),
        legend.position = 'none') +
  stat_anova_test(method = 'one_way')

#----------------------------------------------------
library(survminer)
#Break scores into tertiles (or other ntiles)
q_up = quantile(gsva_h$up, probs = seq(0, 1, by = 1/3), na.rm = TRUE)

#Assign each patient to a tertile
gsva_h$up_bin = cut(gsva_h$up, breaks = q_up, include.lowest = TRUE, labels = FALSE)

#Set the tertiles as factor and establish reference level as the first or the lowest tertile
gsva_h$up_bin = factor(gsva_h$up_bin, levels = c(1, 2, 3))
levels(gsva_h$up_bin)

#Perform regression using the first tertile as the reference level
reg_bin = coxph(s ~ up_bin + age + sex, data = gsva_h)
summary(reg_bin)

#Visualize
fit = survfit(s ~ up_bin, data = gsva_h)
fit #get difference in median survival by computing difference between bins (i.e. bin 3 - bin 2)
p = ggsurvplot(fit, data = gsva_h, risk.table = TRUE, surv.median.line = 'hv', size=0.6)

#Modify the plot by changing strata names and add result from multivariate Cox model
p$plot + 
  labs(x = 'Time (Days)', y = 'Progression Free Probability') +
  scale_color_manual(
    values = c(
      'up_bin=1' = "green3",
      'up_bin=2' = "orange",
      'up_bin=3' = "red"
    ), labels = c("Low Risk", "Medium Risk", "High Risk")) +
  geom_text(aes(x = 1700, y = 0.8, label = 'Median Survival Difference Low-High Risk: 621 Days\nAdjusted HR: 3.3\nP = 7.9e-06'), size = 2.8) 
  







#----------------OS----------------------------------------


## Cox Proportional Hazard Regression
library(survival)

survival = cli$PAAD_survival.txt
survival = survival[match(gsva_h$sampleID, survival$sample), ]
sOS = Surv(time = survival$OS.time, event = survival$OS)

regOS = coxph(sOS ~ up + down, data = gsva_h)
summary(regOS)


## Multivariate Regression
gsva_h = cbind(gsva_h, 
               age = primary$age_at_initial_pathologic_diagnosis, 
               sex = primary$gender,
               grade= primary$neoplasm_histologic_grade)

m_regOS = coxph(sOS ~ up + age + sex + grade, data = gsva_h)
summary(m_regOS)



library(survminer)
#Break scores into tertiles (or other ntiles)
q_up = quantile(gsva_h$up, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
gsva_h$up_bin = cut(gsva_h$up, breaks = q_up, include.lowest = TRUE, labels = FALSE)

gsva_h$up_bin = factor(gsva_h$up_bin, levels = c(1, 2, 3))
levels(gsva_h$up_bin)

#Perform regression using the first tertile as the reference level
reg_binOS = coxph(sOS ~ up_bin + age + sex, data = gsva_h)
summary(reg_binOS)


ggadjustedcurves(
  fit = reg_binOS,
  data = gsva_h,
  variable = "up_bin",
  legend.title = "Risk Group",
  size = 0.6
) +
  labs(x = "Time (Days)", y = "Adjusted Overall Survival") +
  scale_color_manual(
    values = c(
      '1' = "green3",
      '2' = "orange",
      '3' = "red"
    ), labels = c("Low Risk", "Medium Risk", "High Risk")) +
  geom_text(aes(x = 1400, y = 0.8, label = 'Adjusted Low-High Risk HR: 2.7\nP = 1.3e-04'), size = 3,  color = "black")




#Visualize
fitOS = survfit(sOS ~ up_bin, data = gsva_h)
fitOS #get difference in median survival by computing difference between bins (i.e. bin 3 - bin 2)

summary(fitOS, times = c(365, 1095, 1825)) 

pOS = ggsurvplot(fitOS, data = gsva_h, risk.table = TRUE, surv.median.line = 'hv')

#Modify the plot by changing strata names and add result from multivariate Cox model
pOS$plot + 
  labs(x = 'Time (Days)', y = 'Overall Survival') +
  scale_color_manual(
    values = c(
      'up_bin=1' = "green3",
      'up_bin=2' = "orange",
      'up_bin=3' = "red"
    ), labels = c("Low Risk", "Medium Risk", "High Risk")) +
  geom_text(aes(x = 1400, y = 0.8, label = 'Adjusted HR: 3.6\nP = 4.9e-06'), size = 3)


#----------------------GENE PATHWAYS----------------------


hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets1 <- split(hallmark$gene_symbol, hallmark$gs_name)

exp.matrix <- as.matrix(exp.data[,-1])
rownames(exp.matrix) <- exp.data$sample

library(GSVA)

gsva_obj1 = gsvaParam(
  exprData = exp.matrix, 
  geneSets = gene_sets1, 
  #minSize = 5, maxSize = 500, 
  kcdf = 'Gaussian'
)
gsva_h1 = gsva(gsva_obj1) %>% t() %>% as.data.frame()
gsva_h1$sampleID = rownames(gsva_h1)


library(pheatmap)



#pheatmap(as.matrix(gsva_h1),
 #        scale = "row",   # important for expression data
  #       clustering_distance_rows = "euclidean",
   #      clustering_distance_cols = "euclidean",
    #     clustering_method = "complete")

clinical = cli$PAAD_clinicalMatrix
primary1 = subset(clinical, sample_type == 'Primary Tumor')
#Matching sampleIDs between clinical data matrix and gsva result matrix
gsva_h1 <- gsva_h1[match(primary1$sampleID, gsva_h1$sampleID),]
gsva_h1 <- gsva_h1[-which(is.na(gsva_h1$sampleID)),] ## removing mismatched values

primary1 = primary1[match(gsva_h1$sampleID, primary1$sampleID), ] ## removing mismatched values


library(survival)

survival = cli$PAAD_survival.txt
survival = survival[match(gsva_h1$sampleID, survival$sample), ]
s = Surv(time = survival$PFI.time, event = survival$PFI)
reg_path = coxph(s ~ HALLMARK_MYC_TARGETS_V1, data = gsva_h1)
summary(reg_path)


## Multivariate Regression
gsva_h1 = cbind(gsva_h1, 
               age = primary1$age_at_initial_pathologic_diagnosis, 
               sex = primary1$gender)

m_reg1 = coxph(s ~ HALLMARK_MYC_TARGETS_V1 + age + sex, data = gsva_h1)
summary(m_reg1)

library(survminer)
#Break scores into tertiles (or other ntiles)
q_up1 = quantile(gsva_h1$HALLMARK_MYC_TARGETS_V1, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
gsva_h1$up_bin = cut(gsva_h1$HALLMARK_MYC_TARGETS_V1, breaks = q_up1, include.lowest = TRUE, labels = FALSE)

#Set the tertiles as factor and establish reference level as the first or the lowest tertile
gsva_h1$up_bin = factor(gsva_h1$up_bin, levels = c(1, 2, 3))
levels(gsva_h1$up_bin)

#Perform regression using the first tertile as the reference level
reg_bin1 = coxph(s ~ up_bin + age + sex, data = gsva_h1)
summary(reg_bin1)

#Visualize
fit1 = survfit(s ~ up_bin, data = gsva_h1)
fit1 #get difference in median survival by computing difference between bins (i.e. bin 3 - bin 2)
p1 = ggsurvplot(fit1, data = gsva_h1, risk.table = TRUE, surv.median.line = 'hv')

#Modify the plot by changing strata names and add result from multivariate Cox model
p1$plot + 
  labs(x = 'Time (Days)', y = 'Progression Free Probability') +
  scale_color_manual(
    values = c(
      'up_bin=1' = "green3",
      'up_bin=2' = "orange",
      'up_bin=3' = "red"
    ), labels = c("Low Risk", "Medium Risk", "High Risk")) +
  geom_text(aes(x = 1400, y = 0.8, label = 'MYC'), size = 2.8)


analyze_pathway <- function(gsva_h1, pathway_col) {
  
  q_up1 = quantile(gsva_h1[[pathway_col]], probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  
  gsva_h1$up_bin = cut(gsva_h1[[pathway_col]], breaks = q_up1, include.lowest = TRUE, labels = FALSE)
  
  gsva_h1$up_bin = factor(gsva_h1$up_bin, levels = c(1, 2, 3))
  
  reg_bin1 = coxph(s ~ up_bin + age + sex, data = gsva_h1)
  
  fit1 = survfit(s ~ up_bin, data = gsva_h1)
  
  p1 = ggsurvplot(fit1, data = gsva_h1, risk.table = TRUE, surv.median.line = 'hv')
  
  p1$plot + 
    labs(x = 'Time (Days)', y = 'Progression Free Probability') +
    scale_color_manual(
      values = c(
        'up_bin=1' = "green3",
        'up_bin=2' = "orange",
        'up_bin=3' = "red"
      ),
      labels = c("Low Risk", "Medium Risk", "High Risk")
    ) +
    geom_text(aes(x = 1400, y = 0.8, label = pathway_col), size = 2.8)
  
  return(list(
    summary = summary(reg_bin1),
    plot = p1$plot
  ))
}


result <- analyze_pathway(gsva_h1, "HALLMARK_G2M_CHECKPOINT")

result





#----------------------------

analyze_pathway <- function(gsva_h1, pathway_col) {
  
  # Ensure numeric (prevents your earlier error)
  gsva_h1[[pathway_col]] <- as.numeric(as.character(gsva_h1[[pathway_col]]))
  
  # Create tertiles
  q_up1 <- quantile(gsva_h1[[pathway_col]], probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  
  gsva_h1$up_bin <- cut(gsva_h1[[pathway_col]],
                        breaks = q_up1,
                        include.lowest = TRUE,
                        labels = c("Low", "Med", "High"))
  
  # Cox model
  reg_bin1 <- coxph(s ~ up_bin + age + sex, data = gsva_h1)
  
  summ <- summary(reg_bin1)
  
  # 🔥 Extract HR for HIGH vs LOW
  # (row name depends on factor coding → check safely)
  row_name <- grep("High", rownames(summ$coefficients), value = TRUE)
  
  hr <- summ$coefficients[row_name, "exp(coef)"]
  pval <- summ$coefficients[row_name, "Pr(>|z|)"]
  
  return(data.frame(
    pathway = pathway_col,
    HR = hr,
    pvalue = pval
  ))
}

selected_pathways <- gsea %>%
  arrange(NES) %>%   # ascending → most negative first
  slice_head(n = 10)

pathways <- selected_pathways$pathway


results <- do.call(rbind, lapply(pathways, function(pw) {
  analyze_pathway(gsva_h1, pw)
}))

results_sorted <- results[order(-results$HR), ]

head(results_sorted, 10)


write.csv(results_sorted, "pathway_HR_results.csv", row.names = FALSE)

library(knitr)

kable(head(results_sorted, 10))

library(gt)

results_sorted %>%
  slice(1:10) %>%
  gt()