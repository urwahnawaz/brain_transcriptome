## Cleaning and formatting GTEx metadata 

## Loading libraries 
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "purrr")
libsLoaded = lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})




add_feature = function(feature_column, features){
  as.vector(sapply(feature_column, function(x){
    names(features)[sapply(features, function(f) x %in% f)]})) 
}


# Function to preprocess raw bulk data  
process_gtex = function(GTEx.path, outDIR){
  
  source("def_stages.R") 
  
  ## List all files 
  attributes = list.files(GTEx.path, full.names = TRUE, pattern = "\\SampleAttributesDS.txt") # Sample attributes contains sample level information
  phenotype = list.files(GTEx.path, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt") # Phenotype level information related to each donor 
  exp = list.files(GTEx.path, full.names = TRUE, pattern = "\\gene_tpm.gct.gz") # File used for expression matrix 
  
  ## Load all files 
  gtex.md.attributes = read.table(attributes, sep = "\t", header=TRUE, fill=TRUE, quote="")
  gtex.md.pheno = read.table(phenotype, sep = "\t", header=TRUE, fill=TRUE, quote="")
  gtex.annot = read.csv("../annotations/GTEx-metadata-annot.csv", header= TRUE) # Contains custom defined column names for BITHub
  gtex.exp = read.delim(exp, skip = 2)
  
  
  ## Join the phneotype data with the sample information 
  message("To consolidate the sample and phenotype information, we are now joining both SubjectPhenotype and SampleAttribute files")
  message(paste0("Now joining .."))
  
  x = transpose(strsplit(gtex.md.attributes$SAMPID, split="-", fixed=TRUE))
  SampleName=paste(x[[1]], x[[2]], sep="-")
  gtex.md.attributes$SUBJID <- SampleName
  pheno = gtex.md.pheno[match(SampleName, gtex.md.pheno$SUBJID),]
  gtex.md.attributes %<>%   inner_join(pheno, by = "SUBJID") %>% distinct()
  
  message("Merging complete!")
  
  
  message("Changing GTEx column names to BITHub column names")
  ## Changing the GTEx column names to our defined column names for BITHub 
  ## For more information, visit brain_transcriptome/annotations/GTEX-metadata-annot.csv
  colnames(gtex.md.attributes) = gtex.annot$BITColumnName[match(colnames(gtex.md.attributes), gtex.annot$OriginalMetadataColumnName)]
  
  message("Subsetting brain only samples from the metadata file")
  gtex.md.attributes %<>%
    dplyr::filter(grepl("Brain", Structure)) 
  
  message("Adding structure and region annotations")
  
  # structure_acronym and regions lists get generated when running source("def_stages.R")
  
  gtex.md.attributes$StructureAcronym = add_feature(gtex.md.attributes$Structure, structure_acronym) 
  
  gtex.md.attributes$Regions = add_feature(gtex.md.attributes$StructureAcronym, regions)
  
  gtex.md.attributes$Period <- "Postnatal" # All GTEx samples are postnatal 
  gtex.md.attributes$Diagnosis <- "Control" # All GTEx samples have normal brain phenotypes (no neuro-related illnesses reported)
  gtex.md.attributes$AgeInterval <- paste(gtex.md.attributes$AgeInterval, "yrs", sep = "") #Add yrs to match consistency of other datasets
  
  
  # Changing numeric annotations to M or F
  gtex.md.attributes$Sex <- gsub("1", "M", gtex.md.attributes$Sex) 
  gtex.md.attributes$Sex <- gsub("2", "F", gtex.md.attributes$Sex)
  
  
  message("Now pre-processing the expression file")
  colnames(gtex.exp) <- gsub("\\.", "-", colnames(gtex.exp)) # Changing expression file names to match metadata SampleIDs
  gtex.exp %<>% column_to_rownames("Name")
  
  message("Subsetting 'Brain' only samples from the GTEx expression matrix")
  
  
  gtex.brain = gtex.exp %>% 
    dplyr::select(contains(gtex.md.attributes$SampleID))
  message(paste0("Samples subsetted - Exp matrix contains ", ncol(gtex.brain), " samples"))
  
  rownames(gtex.brain) <- sub("\\.[0-9]*$", "", rownames(gtex.brain)) # Removing versions from the Ensembl gene IDs
  
  # Ensuring both samples in metadata and exp file match
  gtex.md.attributes = gtex.md.attributes[which(gtex.md.attributes$SampleID %in% colnames(gtex.brain)),] 
  
  gtex.md.attributes = apply(gtex.md.attributes,2,as.character)
  
  ## Saving all files 
  ## Find outDIR 
  if(file.exists(outDIR)) {
    message("Results directory exists")} else {
      message("Creating results directory")
      dir.create(outDIR, showWarnings = FALSE)}  
  message(paste0("Saving results to ", outDIR))
  write.csv(gtex.md.attributes, paste0(outDIR, "/GTEx-metadata.csv"))
  write.csv(gtex.brain, paste(outDIR, "/GTEx-exp.csv"))
  
  
  
}


gtex = process_gtex(GTEx.path = "../../Data/tmp_data/", outDIR = "../../Data/tmp_data/GTEx/")
gtex.file = read.csv("../../Data/tmp_data/GTEx/GTEx-metadata.csv", header= TRUE)
head(gtex.file)


# BrainSpan 

process_bspan = function(bspan.path, outDIR){
  
  source("def_stages.R")
  columns.metadata = read.csv(file.path(bspan.path, "columns_metadata.csv"), header = TRUE)
  counts.matrix = read.csv(file.path(bspan.path, "expression_matrix.csv"), header= FALSE, row.names= 1)
  rows.metadata = read.csv(file.path(bspan.path, "rows_metadata.csv"))
  md.annot = read.csv()
  bspan.mRIN = readxl::read_xlsx("../BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/mRIN/ncomms8816-s2.xlsx", skip = 1) %>% 
    as.data.frame()
  
  
  columns.metadata$stage = add_feature(columns.metadata$age, stages)
  columns.metadata$stage = factor(columns.metadata$stage, levels(columns.metadata$stage)[match(order.stages, levels(columns.metadata$stage))])
  columns.metadata$Regions = add_feature(columns.metadata$structure_acronym, regions)
  column.metadata$Age_Interval = add_features(columns.metadata$age, age_interval)
  columns.metadata$Age_Interval = factor(columns.metadata$age_interval, levels(columns.metadata$age_interval)[
    match(order.intervals, levels(columns.metadata$age_interval))])
  
  
  ## Conversion to numeric ages 
  message("Now converting to numeric ages")
  columns.metadata$age <-  gsub(" ","_", columns.metadata$age)
  columns.metadata$AgeNumeric[grepl("pcw", columns.metadata$age, ignore.case = TRUE)]<-
    columns.metadata$age[grepl("pcw", columns.metadata$age)] %>%  str_remove("_pcw")%>% 
    as.numeric() %>% `-` (40) %>% divide_by(52)
  columns.metadata$AgeNumeric[grepl("mos", columns.metadata$age, ignore.case = TRUE)] <- 
    columns.metadata$age[grepl("_mos", columns.metadata$age)] %>%  str_remove("_mos") %>% 
    as.numeric() %>% divide_by(12)
  columns.metadata$AgeNumeric[grepl("yrs", columns.metadata$age, ignore.case = TRUE)] <- 
    columns.metadata$age[grepl("_yrs", columns.metadata$age)] %>%  str_remove("_yrs") %>% 
    as.numeric 
  
  
  columns.metadata$Diagnosis <- "Control"
  
  columns.metadata <- columns.metadata %>% 
    dplyr::mutate(SampleID = paste(donor_id, age, structure_acronym, stage, sep = "_"))
  
  columns.metadata %<>% 
    as.data.frame() %>%
    dplyr::select("SampleID", everything())
  
  
  columns.metadata$age_for_mRIN <-  gsub("_","", columns.metadata$age)
  columns.metadata$donor_name <-  gsub("\\.","_", columns.metadata$donor_name)
  columns.metadata = columns.metadata %>% 
    dplyr::mutate(mRIN_id = paste(donor_name, age_for_mRIN,gender, structure_acronym, sep = "//"))
  head(columns.metadata$mRIN_id)
  
  
  bspan.mRIN %<>% dplyr::select("sample name", "mRIN", "z-score", "P-value")
  
  columns.metadata = merge(columns.metadata, bspan.mRIN, by.x="mRIN_id", by.y = "sample name", all = TRUE, no.dups = TRUE)
  columns.metadata = columns.metadata[!is.na(columns.metadata$SampleID),]
  
  columns.metadata = columns.metadata[order(as.numeric(as.character(columns.metadata$column_num))),]
  columns.metadata %<>% dplyr::select(-mRIN_id)
  
  colnames(counts.matrix) <- columns.metadata$SampleID
  rownames(counts.matrix) <- rows.metadata$ensembl_gene_id
  counts.matrix %<>% rownames_to_column("EnsemblID")
}


## progress bar 

temp <- list.files(pattern="*\\.tsv$") 
pb <- progress_bar$new(format = " progress [:bar] :percent eta: :eta", # add
                       total = length(temp), clear = FALSE, width= 60) # add
test_data <- lapply(temp,function(x){
  pb$tick()                                                            # add
  read.csv(file = x,
           sep ="\t",
           fill = TRUE,
           quote='', 
           header = FALSE 
  )[ ,c(287, 288, 289, 290, 291, 292, 293, 304, 370, 661, 662, 812, 813,994, 995, 1002)]
})
