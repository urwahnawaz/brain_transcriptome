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



process_bspan = function(bspan.path, outDIR){
  
  source("def_stages.R")
  
  ## Load all files 
  
  columns.metadata = read.csv(file.path(bspan.path, "columns_metadata.csv"), header = TRUE)
  counts.matrix = read.csv(file.path(bspan.path, "expression_matrix.csv"), header= FALSE, row.names= 1)
  rows.metadata = read.csv(file.path(bspan.path, "rows_metadata.csv"))
  md.annot = read.csv("../../../BrainData/Bulk/FormattedData/BrainSpan/BrainSpan-metadata-annot.csv", header= TRUE)
  bspan.mRIN = readxl::read_xlsx("../annotations/ncomms8816-s2.xlsx", skip = 1) %>% 
    as.data.frame()
  
  message("Changing BrainSpan column names to BITHub column names")
  colnames(columns.metadata) = md.annot$BITColumnName[match(colnames(columns.metadata), md.annot$OriginalMetadataColumnName)]
  
  message("Adding additional metadata information")
  columns.metadata$Stage = add_feature(columns.metadata$Age, stages)
  columns.metadata$Regions = add_feature(columns.metadata$StructureAcroymn, regions)
  columns.metadata$AgeInterval = add_feature(columns.metadata$Age, age_intervals)
  columns.metadata$Diagnosis <- "Control"
  
  
  ## Conversion to numeric ages 
  message("Now ages to numeric ages")
  columns.metadata$Age <-  gsub(" ","_", columns.metadata$Age)
  columns.metadata$AgeNumeric[grepl("pcw", columns.metadata$Age, ignore.case = TRUE)]<-
    columns.metadata$age[grepl("pcw", columns.metadata$Age)] %>%  str_remove("_pcw")%>% 
    as.numeric() %>% `-` (40) %>% divide_by(52)
  columns.metadata$AgeNumeric[grepl("mos", columns.metadata$Age, ignore.case = TRUE)] <- 
    columns.metadata$age[grepl("_mos", columns.metadata$Age)] %>%  str_remove("_mos") %>% 
    as.numeric() %>% divide_by(12)
  columns.metadata$AgeNumeric[grepl("yrs", columns.metadata$Age, ignore.case = TRUE)] <- 
    columns.metadata$age[grepl("_yrs", columns.metadata$Age)] %>%  str_remove("_yrs") %>% 
    as.numeric 
  
  message("Creating SampleIDs using Donor ID, Age, Structure Acronym and Stage ")
  columns.metadata <- columns.metadata %>% 
    dplyr::mutate(SampleID = paste(DonorID, Age, StructureAcromyn, Stage, sep = "_"))
  
  columns.metadata %<>% 
    as.data.frame() %>%
    dplyr::select("SampleID", everything())
  
  message("Adding mRIN information - based from Feng et al (2015)")
  columns.metadata$age_for_mRIN <-  gsub("_","", columns.metadata$Age)
  columns.metadata$DonorName <-  gsub("\\.","_", columns.metadata$DonorName)
  columns.metadata = columns.metadata %>% 
    dplyr::mutate(mRIN_id = paste(DonorName, age_for_mRIN,Sex ,StructureAcroymn, sep = "//"))
  
  
  bspan.mRIN %<>% dplyr::select("sample name", "mRIN", "z-score", "P-value")
  
  columns.metadata = merge(columns.metadata, bspan.mRIN, by.x="mRIN_id", by.y = "sample name", all = TRUE, no.dups = TRUE)
  columns.metadata = columns.metadata[!is.na(columns.metadata$SampleID),]
  
  columns.metadata = columns.metadata[order(as.numeric(as.character(columns.metadata$column_num))),]
  columns.metadata %<>% dplyr::select(-mRIN_id, -age_for_mRIN)
  
  message("Preprocessing of metadata file complete")
  message("Now processing expression file")
  
  colnames(counts.matrix) <- columns.metadata$SampleID
  rownames(counts.matrix) <- rows.metadata$ensembl_gene_id
  counts.matrix %<>% rownames_to_column("EnsemblID")
  
  
  columns.metadata = apply(columns.metadata,2,as.character)
  
  ## Saving all files 
  ## Find outDIR 
  if(file.exists(outDIR)) {
    message("Results directory exists")} else {
      message("Creating results directory")
      dir.create(outDIR, showWarnings = FALSE)}  
  message(paste0("Saving results to ", outDIR))
  write.csv(columns.metadata, paste0(outDIR, "/BrainSpan-metadata.csv"))
  write.csv(counts.matrix, paste(outDIR, "/BrainSpan-exp.csv"))
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
