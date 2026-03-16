#Set-up
    library("here")
    source(here("src/setup/header.R"))
    count_files = list.files("data/raw/counts/", full.names = TRUE)
    annotatedIDs = list.files("data/raw/metadata/", full.names = TRUE)

#Data cleaning
    dat = sapply(
      count_files, 
      function(x){
        read.delim(x, header = F) 
        }
    )
    dat = do.call('cbind', dat)               #bind dataframe
    rownames(dat) <- dat[,1]                  #set rownames
    dat <- dat[, dat[1, ] != "Geneid"]        #remove redundant columns
    names(dat) <- substr(names(dat), 18, 20)  #set names
    dat<- dat[-1,]                            #remove sample labels
    dat<- dat[,-18]                           #remove gene lengths
    dat[] = lapply(dat, as.integer)           #convert to numeric
    dat$KP5 <- dat$KP5 + dat$KP5.1            #average technical replicates
    dat <- dat[,!colnames(dat) %in% "KP5.1"]  #remove extra col

#Make Annotations
    #geneIDs = read.delim(paste0(annotatedIDs), header = T)
    cond_code <- substr(colnames(dat), 1, 1) #Condition from 1st letter
    condition <- dplyr::case_when(
      cond_code == "C" ~ "Control",
      cond_code == "K" ~ "Treatment",
      TRUE ~ NA_character_
    )
    cell_line <- substr(colnames(dat), 2, 3) #Cell lines from 2-3 character
    annot <- data.frame(
      condition   = factor(condition, levels = c("Control", "Treatment")),
      cell_line   = factor(cell_line),
      stringsAsFactors = FALSE
    )
    rownames(annot) <- colnames(dat)         #Add rowwnames

#Write Files
    out_dir <- file.path("data/processed/")
    # write processed data
    write.table(
      dat,
      file.path(out_dir, "counts.tsv"),
      sep = "\t",
      col.names = TRUE,
      row.names = TRUE,
      quote = FALSE
    )
    # write sample annotation
    write.table(
      annot,
      file.path(out_dir, "sample_annotations.txt"),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE,
    )
    message(
      "Complete, outputs written to: ",
      out_dir,
      "\n"
    )
    rm(list = ls())
