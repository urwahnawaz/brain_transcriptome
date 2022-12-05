## Cleaning and formatting GTEx metadata 

## Loading libraries 
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "purrr")
libsLoaded = lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})




add_feature = function(feature_column, features){
  as.vector(sapply(feature_column, function(x){
  names(features)[sapply(features, function(f) x %in% f)]})) 
}


# Function to preprocess GTEx data 
process_gtex = function(GTEx.path, outDIR){
  
  source("def_stages.R")
  ## Needs a mechanism to list all files and make sure all are included 
  # match patterns 
  attributes = list.files(GTEx.path, full.names = TRUE, pattern = "\\SampleAttributesDS.txt")
  phenotype = list.files(GTEx.path, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt")
  #exp = list.file(GTEx_path, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt")
  ## Load all files 
  gtex.md.attributes = read.table(attributes, sep = "\t", header=TRUE, fill=TRUE, quote="")
  gtex.md.pheno = read.table(phenotype, sep = "\t", header=TRUE, fill=TRUE, quote="")
  gtex.annot = read.csv("../annotation_files /GTEx-metadata-annot.csv", header= TRUE)
  #gtex.exp = read.delim()

  
  ## Join the phneotype data with the sample information 
  message("To have a complete picture of the GTEx metadata, we are now joining both SubjectPhenotype and SampleAttribute files")
  message(paste0("Now joining .."))
  x = transpose(strsplit(gtex.md.attributes$SAMPID, split="-", fixed=TRUE))
  SampleName=paste(x[[1]], x[[2]], sep="-")
  gtex.md.attributes$SUBJID <- SampleName
  pheno = gtex.md.pheno[match(SampleName, gtex.md.pheno$SUBJID),]
  gtex.md.attributes %<>%   inner_join(pheno, by = "SUBJID") %>% distinct()
  
  message("Complete!")
  

  message("Changing GTEx column names to BITHub column names")
  colnames(gtex.md.attributes) = gtex.annot$BITColumnName[match(colnames(gtex.md.attributes), gtex.annot$OriginalMetadataColumnName)]
  
  
  gtex.md.attributes %<>% 
    dplyr::filter(grepl("Brain", Structure)) 
  
  message("Adding structure and region annotations")
  
  gtex.md.attributes$StructureAcronym = add_feature(gtex.md.attributes$Structure, structure_acronym)
  
  gtex.md.attributes$Regions = add_feature(gtex.md.attributes$StructureAcronym, regions)
  
  gtex.md.attributes$Period <- "Postnatal"
  gtex.md.attributes$Diagnosis <- "Control"
  gtex.md.attributes$AgeInterval <- paste(gtex.md.attributes$AgeInterval, "yrs", sep = "")
  
  
  gtex.md.attributes$Sex <- gsub("1", "M", gtex.md.attributes$Sex)
  gtex.md.attributes$Sex <- gsub("2", "F", gtex.md.attributes$Sex)
  
  ## Find outDIR 
  if(file.exists(outDIR)) {
    message("Results directory exists")} else {
      message("Creating results directory")
      dir.create(outDIR, showWarnings = FALSE)}  
  message(paste0("Saving results to ", outDIR))
  write.csv(gtex.md.attributes, paste0(outDIR, "GTEx-metadata.csv"))
  
  

}


gtex = process_gtex(GTEx.path = "../../Data/tmp_data/")
gtex

## testing code 
GTEx.path = file.path("../../Data/tmp_data")
attributes = list.files(GTEx.path, full.names = TRUE, pattern = "\\SampleAttributesDS.txt")
phenotype = list.files(GTEx.path, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt")
attributes

# Test the code here 
## Find outDIR 
#  if(file.exists(outDIR)) {
#   message("Results directory exists")} else {
#      message("Creating results directory")
#     dir.create(outDIR, showWarnings = FALSE)}  
#  message(paste0("Saving results to ", outDIR))
#  write.csv(gtex.md.attributes, paste0(outDIR, "GTEx-metadata.csv"))


gtex.md.attributes = read.table("../../Data/tmp_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
                                sep = "\t", header=TRUE, fill=TRUE, quote="")
gtex.md.pheno = read.table("../../Data/tmp_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", 
                           sep="\t", header=TRUE, fill=TRUE, quote="")

gtex.annot = read.csv("../../Data/FormattedData/Gtex/GTEx-metadata-annot.csv", header = TRUE)
structures = read.csv("../annotation_files /structure_acronyms.csv", header= TRUE)


## Combine both data frames together 
x<- transpose(strsplit(gtex.md.attributes$SAMPID, split="-", fixed=TRUE))
SampleName=paste(x[[1]], x[[2]], sep="-")

gtex.md.attributes$SUBJID <- SampleName
pheno <- gtex.md.pheno[match(SampleName, gtex.md.pheno$SUBJID),]
pheno
gtex.md.attributes %<>% 
  inner_join(pheno, by = "SUBJID") %>% distinct()
gtex.md.attributes


## Add gtex annotations 
colnames(gtex.md.attributes) = gtex.annot$BITColumnName[match(colnames(gtex.md.attributes), gtex.annot$OriginalMetadataColumnName)]

gtex.md.attributes %<>% 
  dplyr::filter(grepl("Brain", Structure)) 

gtex.md.attributes$StructureAcronym = as.vector(sapply(gtex.md.attributes$Structure, function(x){
  names(structure_acronym)[sapply(structure_acronym, function(structure) x %in% structure)]})) 

gtex.md.attributes$StructureAcronym = add_feature(gtex.md.attributes$Structure, structure_acronym)

## Add other annotations 
gtex.md.attributes$StructureAcronym = structures$BITHub[match(gtex.md.attributes$Structure, structures$OriginalName)]

gtex.md.attributes$Regions = as.vector(sapply(gtex.md.attributes$StructureAcronym, function(x){
  as.vector(names(regions))[sapply(regions, function(region) x %in% region)]}))

gtex.md.attributes

gtex.md.attributes$Period <- "Postnatal"
gtex.md.attributes$Diagnosis <- "Control"
gtex.md.attributes$AgeInterval <- paste(gtex.md.attributes$AgeInterval, "yrs", sep = "")


gtex.md.attributes$Sex <- gsub("1", "M", gtex.md.attributes$Sex)
gtex.md.attributes$Sex <- gsub("2", "F", gtex.md.attributes$Sex)

head(gtex.md.attributes)