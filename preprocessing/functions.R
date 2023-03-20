libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "purrr")
libsLoaded = lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


source("def_stages.R")

add_feature = function(feature_column, features){
  as.vector(sapply(feature_column, function(x){
    names(features)[sapply(features, function(f) x %in% f)]})) 
}

num_to_round = function(age){
  if (is.na(age)) {
    NaN
  } else if (age >= 2) {
    paste0(round(age), " yrs")
  } else if (age < 0) {
    paste0(round(age * 52 + 40), " pcw")
  } else if (age >= 0 & age < 2) {
    paste0(round(age * 12), " mos")
  }
}



clean_and_format = function(dir ,dataset, outdir){
  
  message(paste0("Now formatting the ", dataset, " dataset"))
  annot= read.csv(paste0("../annotations/", dataset,"-metadata-annot.csv"))
  
  if (dataset == c("GTEx")) {
    message("Retriving files for GTEx")
    
    attributes = list.files(dir, full.names = TRUE, pattern = "\\SampleAttributesDS.txt") # Sample attributes contains sample level information
    phenotype = list.files(dir, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt") # Phenotype level information related to each donor 
    exp = list.files(dir, full.names = TRUE, pattern = "\\gene_tpm.gct.gz") # File used for expression matrix 
    
    md = read_tsv(attributes, col_types = c('.default' = 'c')) %>% 
      filter(SMTS == 'Brain') %>% 
      mutate(SUBJID = sapply(str_split(SAMPID, pattern = "-"), function(x) paste(x[1:2], collapse = '-'))) %>%
      left_join(read_tsv(phenotype, col_types = c('.default' = 'c')))  %>% as.data.frame()
    colnames(md) = annot$BITColumnName[match(colnames(md), annot$OriginalMetadataColumnName)]
    
    md = md %>% mutate(StructureAcronym = add_feature(.$Structure, structure_acronym)) %>% 
      mutate(Regions = add_feature(.$StructureAcronym, regions), 
             AgeInterval = paste(.$Age, "yrs", sep = ""), 
             Diagnosis = "Control", 
             Sex = ifelse(Sex == 1, "M", "F"), 
             Period = "Postnatal") %>% as.data.frame()
    
    exp = read.delim(exp, skip = 2)
    colnames(exp) <- gsub("\\.", "-", colnames(exp)) # Changing expression file names to match metadata SampleIDs
    exp %<>% column_to_rownames("Name")
    
    exp = exp %>% 
      dplyr::select(contains(md$SampleID))
    message(paste0("Samples subsetted - Exp matrix contains ", ncol(exp), " samples"))
    rownames(exp) <- sub("\\.[0-9]*$", "", rownames(exp))
    
    md= md[which(md$SampleID %in% colnames(exp)),] 
    
    
    
    
  } else if (dataset == "BrainSpan"){
    # Loading all files 
    
    columns.metadata = read.csv(file.path(dir, "columns_metadata.csv"), header = TRUE)
    exp = read.csv(file.path(dir, "expression_matrix.csv"), header= FALSE, row.names= 1)
    rows.metadata = read.csv(file.path(dir, "rows_metadata.csv"))
    
    bspan.mRIN = readxl::read_xlsx("../annotations/ncomms8816-s2.xlsx", skip = 1) %>% 
      as.data.frame() %>% dplyr::select("sample name", "mRIN", "z-score", "P-value")
    
    message("Changing BrainSpan column names to BITHub column names")
    colnames(columns.metadata) = annot$BITColumnName[match(colnames(columns.metadata), annot$OriginalMetadataColumnName)]
    
    message("Adding additional metadata information")
  
    md = columns.metadata %>% 
      mutate(Stage = add_feature(.$Age, stages), 
             Regions = add_feature(.$StructureAcronym, regions), 
             AgeInterval = add_feature(.$Age, age_intervals), 
             Diagnosis = "Control", 
             Age = gsub(" ","_", .$Age)) %>%  
      mutate(SampleID = paste(DonorID, Age, StructureAcronym, Stage, sep = "_"), 
             age_for_mRIN = gsub("_", "", .$Age), 
             DonorName = gsub("\\.","_", .$DonorName)) %>% 
      mutate('sample name' = paste(DonorName, age_for_mRIN,Sex ,StructureAcronym, sep = "//")) %>% 
      left_join(bspan.mRIN) %>% as.data.frame() %>% 
      dplyr::select("SampleID", everything())
    
    ## Adding Age Numeric 
    md$AgeNumeric[grepl("pcw", md$Age, ignore.case = TRUE)]<-
      md$Age[grepl("pcw", md$Age)] %>%  str_remove("_pcw")%>% 
      as.numeric() %>% `-` (40) %>% divide_by(52)
    md$AgeNumeric[grepl("mos", md$Age, ignore.case = TRUE)] <- 
      md$Age[grepl("_mos", md$Age)] %>%  str_remove("_mos") %>% 
      as.numeric() %>% divide_by(12)
    md$AgeNumeric[grepl("yrs", md$Age, ignore.case = TRUE)] <- 
      md$Age[grepl("_yrs", md$Age)] %>%  str_remove("_yrs") %>% 
      as.numeric 
    
    md %<>% mutate(Period = ifelse(.$AgeNumeric >= 0, "Postnatal", "Prenatal"), 
                   colname = paste(DonorName, Age, StructureAcronym, sep = "_"))
    
    
    ## Add from additional excel file retrived from Allen Brain Atlas  
    md.excel = read_excel("../annotations/BrainSpan-additional.xlsx",sheet =2, col_names = TRUE, skip =1) %>% 
      as.data.frame() %>% mutate_at(.vars = "AllenInstituteID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>%
      mutate_at(.vars = "Age", .funs = gsub, pattern = "PCW", replacement = "_pcw") %>% 
      mutate_at(.vars = "Age", .funs = gsub, pattern = "M", replacement = "_mos") %>%
      mutate_at(.vars = "Age", .funs = gsub, pattern = "Y", replacement = "_yrs") %>%
      mutate_at(.vars = "Region/Area", .funs = gsub, pattern = "\\/", replacement = "-") %>% 
      mutate(colname = paste(AllenInstituteID, Age, `Region/Area`, sep = "_")) %>% 
      dplyr::select(-c(Agerange, Age, Description))
    
    ## Join the two files 
    md = md %>% 
      left_join(.,md.excel, by = "colname", keep = TRUE)
    
    
    ## Process sheet 2
    md.excel = read_excel("../annotations/BrainSpan-additional.xlsx", sheet = 1, col_names = TRUE) %>% 
      as.data.frame() %>% 
      mutate_at(.vars="Internal ID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>% 
      dplyr::rename("Braincode" = "External ID")
    
    md = merge(md, md.excel, by ="Braincode") 
    md %<>% 
      dplyr::select(-c("Age.y", "colname.y", "colname.x", "Gender", AllenInstituteID,
                     "Region/Area", "age_for_mRIN", "sample name", "Internal ID")) %>% 
      dplyr::rename("Ethnicity"="Ethn.")  
    
    
    md = md[!duplicated(md[,c('column_num')]),]
    md %>% 
      dplyr::arrange(column_num)
    
    
    message("Now adding Sample Information on expression matrix")
    colnames(exp) <- md$SampleID
    rownames(exp) <- rows.metadata$ensembl_gene_id
    exp %<>% rownames_to_column("EnsemblID")
    
    
  
    
    
    
  } else if (dataset == "BrainSeq"){
    load(file.path(dir, "rse_gene_unfiltered.Rdata"), envir = .GlobalEnv)
    load(file.path(dir,"methprop_pd.Rdata"), envir = .GlobalEnv)
    x = rse_gene@colData 
    x <- as.data.frame(x)
    x <- as.data.frame(t(x))
    replicated <- colnames(x)[grep(",", x["SAMPLE_ID",])]
    y <- as.list(x)
    y[replicated] <- lapply(y[replicated], function(z) {
      # which variables to merge
      to.weight <- which(sapply(z, length) > 1 & sapply(z, class) %in% c("numeric", "integer"))
      
      # weighting of the merge
      weighting <- z$numReads # total reads
      weighting <- weighting / sum(weighting) # rather than a straight average, it's based on the number of reads
      
      # apply weighting
      z[to.weight] <- lapply(z[to.weight], function(zz) {
        if (length(weighting) == length(zz)) {
          sum(weighting * zz)
        } else {
          NaN
        }
        
      })
      
      # quickly fix character variables
      char <- which(sapply(z, length) > 1 & sapply(z, class) == "character")
      z[char] <- lapply(z[char], function(zz) {
        paste(zz, collapse = " & ")
      })
      
      return(z)
    })
    
    w <- lapply(y, as.data.frame)
    w <- do.call("rbind", w)
    
    comp <- as.data.frame(pd)
    comp <- comp[,57:64]
    m <- match(rownames(comp), rownames(w)) # they are
    md<- cbind(w, comp[m,])
    colnames(md) = annot$BITColumnName[match(colnames(md), annot$OriginalMetadataColumnName)]
    
    
    # Adding features 
    md %<>% mutate(Period = ifelse(.$Age > 0, "Postnatal", "Prenatal"), 
                   StructureAcronym = gsub("HIPPO", "HIP", .$StructureAcronym),
                   Diagnosis = gsub("Schizo", "Schizophrenia", .$Diagnosis)) %>%
      mutate(Regions = add_feature(.$StructureAcronym, regions)) %>% 
      mutate(Age_rounded = as.character(sapply(na.omit(.$AgeNumeric), num_to_round))) %>% as.data.frame() %>%
      mutate(AgeInterval = as.character(add_feature(.$Age_rounded, age_intervals))) %>% 
      dplyr::select(-Age_rounded) %>%
      dplyr::select("SampleID", everything()) %>%
      as.data.frame()
    
    exp = rse_gene@assays@.xData$data$rpkm
    rownames(exp) <- sub("\\.[0-9]*$", "", rownames(exp))
    
  } else if (dataset == "PsychEncode"){
    
    
    exp = list.files(dir, full.names = TRUE, pattern = "\\Gene_expression_matrix_TPM.txt") %>% 
      read.table(., header=TRUE, row.names = 1, check.names = FALSE)
    md = list.files(dir, full.names = TRUE, pattern = "Job*") %>% read.csv(., header=TRUE)
    comp = list.files(dir, full.names = TRUE, pattern = "\\Cell_fractions*") %>% read_excel() %>%
      as.data.frame() %>% 
      column_to_rownames("CellType")
    
    colnames(md) = annot$BITColumnName[match(colnames(md), annot$OriginalColumnName)]
    # Fix existing columns 
    comp = comp[,-1] 
    comp = as.data.frame(t(comp))
    m <- match(md$SampleID, rownames(comp))
    md <- cbind(md, comp[m,])
    
    
    # PCW to age numeric 
    md$AgeNumeric[grepl("PCW", md$AgeNumeric, ignore.case = TRUE)] = md$AgeNumeric[grepl("PCW", md$AgeNumeric)] %>%  
      str_remove("PCW")%>% 
      as.numeric() %>% `-` (40) %>% divide_by(52)
    
    md$AgeNumeric = gsub("90+", "91", md$AgeNumeric)
    md$AgeNumeric =gsub("\\+", "", md$AgeNumeric)
    
    
    md$AgeNumeric <- as.numeric(as.character(md$AgeNumeric))
    
    md %<>%
      filter(Diagnosis == "Affective Disorder" |
               Diagnosis == "Autism Spectrum Disorder" | 
               Diagnosis == "Bipolar Disorder" |
               Diagnosis == "Control" |
               Diagnosis == "Schizophrenia") %>% 
      mutate(Structure = c("Dorsolateral Prefrontal Cortex"),  ## Adding name of structure
             StructureAcronym = c("DLPFC")) %>%  
      mutate(Period = ifelse(.$AgeNumeric >= 0, "Postnatal", "Prenatal"))  %>%
      mutate(Age_rounded = as.character(sapply(.$AgeNumeric, num_to_round))) %>% as.data.frame() %>%
      mutate(AgeInterval = as.character(add_feature(.$Age_rounded, age_intervals))) %>% 
      mutate(Regions = c("Cortex")) %>% 
      mutate(DonorID = as.character(.$SampleID)) %>%
      dplyr::select(-Age_rounded) %>%
      as.data.frame()
    
    
    exp %<>% rownames_to_column("EnsemblID")
    
  }
  
  
  else {
    message("Invalid dataset")
  }
  
  
  md = apply(md,2,as.character)
  if(file.exists(outdir)) {
    message("Results directory exists")} else {
      message("Creating results directory")
      dir.create(outdir, showWarnings = FALSE)}  
  message(paste0("Saving results to ", outdir))
  write.csv(md, paste0(outdir, "/", dataset, "-metadata.csv"))
  write.csv(exp, paste0(outdir, "/",dataset, "-exp.csv"))
  
  
}



