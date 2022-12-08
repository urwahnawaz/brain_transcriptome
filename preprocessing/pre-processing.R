
## Cleaning and formatting GTEx metadata 
source(c("functions.R", "def_stages.R"))






# Function to preprocess raw bulk data  
process_gtex = function(GTEx.path, outDIR){
  

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
  md.annot = read.csv("../annotations/BrainSpan-metadata-annot.csv", header= TRUE)
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


process_bseq = function(path, outdir){
  source("def_stages.R")
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
  final <- cbind(w, comp[m,])
  
  # Adding features 
  final$Period <- ifelse(final$Age > 0, "Postnatal", "Prenatal")
  final$Region <- gsub("HIPPO", "HIP", final$Region)
  final$Regions = add_feature(final$Region, regions)
  
  
  # Age Intervals 
  final %<>% mutate(AgeInterval =
                      case_when(
                        between(Age, 0, 1.5) ~ "0-5mos", 
                        between(Age,2,5.99) ~ "19mos-5yrs", 
                        between(Age,6,11.99) ~ "6-11yrs", 
                        between(Age, 12,19.99) ~ "12-19yrs",
                        between(Age, 20,29.99) ~ " 20-29yrs",
                        between(Age, 30,39.99) ~ "30-39yrs",
                        between(Age, 40, 49.99) ~ "40-49yrs", 
                        between(Age, 50, 59.99) ~ "50-59yrs", 
                        between(Age, 60, 69.99) ~ "60-69yrs", 
                        between(Age, 70, 79.99) ~ "70-79yrs", 
                        between(Age, 80, 89.99) ~ "80-89yrs", 
                        between(Age,90,99.99) ~ "90-99yrs"))
  
  
  ### Finding intervals manually 
  final$AgeInterval[final$Age <= 0 & final$Age >= 1.5] = "0-5mos"
  final$AgeInterval[final$Age <= 2 & final$Age >= 5.99] = "19mos-5yrs"
  final$AgeInterval[final$Age == -0.5945210] <- "8-9pcw"
  final$AgeInterval[final$Age >= -0.52 & final$Age <= -0.47] <- "13-15pcw"
  final$AgeInterval[final$Age >= -0.47 & final$Age <= -0.42] <- "16-18pcw"
  final$AgeInterval[final$Age >= -0.41 & final$Age <= -0.33] <- "19-24pcw"
  final$AgeInterval[final$Age >= -0.27 & final$Age <= -0.090] <- "25-38pcw"
  
  ## Exp file 
  
  exp = rse_gene@assays@.xData$data$rpkm
  rownames(exp) = sub("\\.[0-9]*$", "", rownames(exp)) 
  
  if(file.exists(outdir)) {
    message("Results directory exists")} else {
      message("Creating results directory")
      dir.create(outdir, showWarnings = FALSE)}  
  message(paste0("Saving results to ", outdir))
  write.csv(final, paste0(outdir, "/BrainSeq-metadata.csv"))
  write.csv(exp, paste0(outdir, "/BrainSeq-exp.csv"))
  
  
}



## progress bar 

#temp <- list.files(pattern="*\\.tsv$") 
#pb <- progress_bar$new(format = " progress [:bar] :percent eta: :eta", # add
 #                      total = length(temp), clear = FALSE, width= 60) # add
#test_data <- lapply(temp,function(x){
#  pb$tick()                                                            # add
#  read.csv(file = x,
#           sep ="\t",
#           fill = TRUE,
#           quote='', 
#           header = FALSE 
#  )[ ,c(287, 288, 289, 290, 291, 292, 293, 304, 370, 661, 662, 812, 813,994, 995, 1002)]
#})
