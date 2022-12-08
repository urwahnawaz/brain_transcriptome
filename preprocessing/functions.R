libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "purrr")
libsLoaded = lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


source("def_stages.R")

add_feature = function(feature_column, features){
  as.vector(sapply(feature_column, function(x){
    names(features)[sapply(features, function(f) x %in% f)]})) 
}


clean_and_format = function(dir ,dataset){
  
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
      left_join(read_tsv(phenotype, col_types = c('.default' = 'c'))) %>%
      mutate(StructureAcronym = add_feature(.$SMTSD, structure_acronym)) %>% 
      mutate(Regions = add_feature(.$StructureAcronym, regions)) %>%  
      mutate(AGE = paste(.$AGE, "yrs", sep = "") ) %>% as.data.frame() 
    colnames(gtex.md) = annot$BITColumnName[match(colnames(gtex.md), annot$OriginalMetadataColumnName)]
    
    md = md %>% 
      mutate(Diagnosis = "Control") %>% 
      mutate(SEX = ifelse(SEX == 1, "M", "F")) %>% 
      mutate(Period = "Postnatal") %>% as.data.frame()
     
    colnames(exp) <- gsub("\\.", "-", colnames(exp)) # Changing expression file names to match metadata SampleIDs
    exp %<>% column_to_rownames("Name")
    
    exp = exp %>% 
      dplyr::select(contains(md$SampleID))
    message(paste0("Samples subsetted - Exp matrix contains ", ncol(exp), " samples"))
    rownames(exp) <- sub("\\.[0-9]*$", "", rownames(exp))
    
    md= md[which(md$SampleID %in% colnames(exp)),] 
    
    
    
    
  } else if (dataset == "BrainSpan"){
    # Loading all files 
    
    columns.metadata = read.csv(file.path(bspan.path, "columns_metadata.csv"), header = TRUE)
    counts.matrix = read.csv(file.path(bspan.path, "expression_matrix.csv"), header= FALSE, row.names= 1)
    rows.metadata = read.csv(file.path(bspan.path, "rows_metadata.csv"))
    
    bspan.mRIN = readxl::read_xlsx("../annotations/ncomms8816-s2.xlsx", skip = 1) %>% 
      as.data.frame() %>% dplyr::select("sample name", "mRIN", "z-score", "P-value")
    
    message("Changing BrainSpan column names to BITHub column names")
    colnames(columns.metadata) = annot$BITColumnName[match(colnames(columns.metadata), annot$OriginalMetadataColumnName)]
    
    message("Adding additional metadata information")
    md = columns.metadata %>% 
      mutate(Stage = add_feature(.$Age, stages), 
             Regions = add_feature(.$StructureAcronym, regions), 
             AgeIntervals = add_feature(.$Age, age_intervals), 
             Diagnosis = "Control") %>% 
      mutate(SampleID = paste(DonorID, Age, StructureAcromyn, Stage, sep = "_"), 
             age_for_mRIN = gsub("_", "", .$Age), 
             DonorName = gsub("\\.","_", .$DonorName)) %>% 
      mutate('sample name' = paste(DonorName, age_for_mRIN,Sex ,StructureAcroymn, sep = "//")) %>% 
      left_join(bspan.mRIN) %>% as.data.frame() %>% 
      dplyr::select("SampleID", everything())
    
    ## Adding Age Numeric 
    md$AgeNumeric[grepl("pcw", md$Age, ignore.case = TRUE)]<-
      md$age[grepl("pcw", md$Age)] %>%  str_remove("_pcw")%>% 
      as.numeric() %>% `-` (40) %>% divide_by(52)
    md$AgeNumeric[grepl("mos", md$Age, ignore.case = TRUE)] <- 
      md$age[grepl("_mos", md$Age)] %>%  str_remove("_mos") %>% 
      as.numeric() %>% divide_by(12)
    md$AgeNumeric[grepl("yrs", md$Age, ignore.case = TRUE)] <- 
      md$age[grepl("_yrs", md$Age)] %>%  str_remove("_yrs") %>% 
      as.numeric 
    
    colnames(exp) <- md$SampleID
    rownames(exp) <- rows.metadata$ensembl_gene_id
    exp %<>% rownames_to_column("EnsemblID")
    
    
    md = apply(md,2,as.character)
      
      
    
  } else if (dataset == "BrainSeq"){
    load(file.path(path, "rse_gene_unfiltered.Rdata"))
    load(file.path(path,"methprop_pd.Rdata"))
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
    
    # Adding features 
    md = md %>% mutate(Period = ifelse(.$Age > 0, "Postnatal", "Prenatal"), 
                             Region = gsub("HIPPO", "HIP", .$Region)) %>%
      mutate(Regions = add_feature(.$Region, regions)) %>% 
      mutate(age_interval = as.character(cut(Age, seq(-1, 100, by = 10)))) %>%
      mutate(AgeInterval = sapply(age_interval, function(i) {
        paste0( as.numeric(gsub("^\\(([-0-9]+),.+", "\\1", i)) + 1,
          "-", as.numeric(gsub(".+,([0-9]+)\\]$", "\\1", i)), "yrs")})) %>% 
      dplyr::select(-age_interval)
    
    
    md$AgeInterval[md$Age >= -0.76 & md$Age <= -0.70] = "1-3pcw"
    md$AgeInterval[md$Age >= -0.701 & md$Age <= -0.62] = "4-7pcw"
    md$AgeInterval[md$Age >= -0.621 & md$Age <= -0.58] = "8-9pcw"
    md$AgeInterval[md$Age >= -0.57 & md$Age <= -0.53] = "10-12pcw"
    md$AgeInterval[md$Age >= -0.52 & md$Age <= -0.48] = "13-15pcw"
    md$AgeInterval[md$Age >= -0.47 & md$Age <= -0.42] = "16-18pcw"
    md$AgeInterval[md$Age >= -0.41 & md$Age <= -0.30] = "19-24pcw"
    md$AgeInterval[md$Age >= -0.29 & md$Age <= -0.038] = "25-38pcw"
    md$AgeInterval[md$Age >= -0.019 & md$Age < 0] = "39-40pcw"
    md$AgeInterval[md$Age >= 0 & md$Age <= 0.49] <- "0-5mos"
    md$AgeInterval[md$Age >= 0.50 & md$Age <= 1.58] <- "6-18mos"
    md$AgeInterval[md$Age >= 1.5833 & md$Age <= 5.99] <- "19mos-5yrs"
    md$AgeInterval[md$Age >= 6 & md$Age <= 11.99] <- "6-11yrs"  
    
    exp = rse_gene@assays@.xData$data$rpkm
    rownames(exp) <- sub("\\.[0-9]*$", "", rownames(exp))
    
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

