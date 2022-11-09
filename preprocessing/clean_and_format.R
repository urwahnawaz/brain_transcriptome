# Loading all libraries
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
"data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

# Path to all files 
path.data <- file.path("../BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/") 
out.path <- file.path("")


columns.metadata <- read.csv(file.path(path.data, "columns_metadata.csv"), header = TRUE)
counts.matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)
rows.metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))

## Age 
order.age <- c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw",
               "4 mos","10 mos",
               "1 yrs","2 yrs","3 yrs","4 yrs","8 yrs",
               "11 yrs","13 yrs","15 yrs","18 yrs","19 yrs",
               "21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")


#### Add stage column
columns.metadata$stage <- as.factor(sapply(columns.metadata$age, function(x) {names(stages)[sapply(stages, function(stage) x %in% stage)]}))
#### Sort factor levels of "stage"
columns.metadata$stage <- factor(columns.metadata$stage, levels(columns.metadata$stage)[match(order.stages, levels(columns.metadata$stage))])
table(columns.metadata$stage)


## Add age interval column 
columns.metadata$age_interval = as.factor(sapply(columns.metadata$age, function(x) {
    names(age_intervals)[sapply(age_intervals,function(age_interval) x %in% age_interval)]}))
columns.metadata$age_interval <-
    factor(columns.metadata$age_interval, levels(columns.metadata$age_interval)[match(order.intervals, levels(columns.metadata$age_interval))])




#### Add period 
#columns.metadata$Period = as.factor(sapply(columns.metadata$age_interval, function(x) {
#    names(period)[sapply(period, function(Period) x %in% period)]}))
#columns.metadata$Period =
#    factor(columns.metadata$Period, levels(columns.metadata$Period)[match(order.period, levels(columns.metadata$Period))])
#table(columns.metadata$Period)
#columns.metadata

columns.metadata$Period[columns.metadata$stage == "s2a"|
                            columns.metadata$stage == "s2b"|
                            columns.metadata$stage == "s3a"|
                            columns.metadata$stage == "s3b"|
                            columns.metadata$stage == "s4"|
                            columns.metadata$stage == "s5"] <- "Prenatal"



columns.metadata$Period[columns.metadata$stage == "s6"|
                            columns.metadata$stage == "s7"|
                            columns.metadata$stage == "s8"|
                            columns.metadata$stage == "s9"|
                            columns.metadata$stage == "s10"|
                            columns.metadata$stage == "s11" |
                            columns.metadata$stage == "s12" |
                            columns.metadata$stage == "s13"] <- "Postnatal"


### Adding regions 
##columns.metadata$Regions = sapply(columns.metadata$structure_acronym, function(x) {
  #  names(regions)[sapply(regions,function(Regions) x %in% regions)]})
#columns.metadata

columns.metadata$Regions[columns.metadata$structure_acronym== "AMY" |
                             columns.metadata$structure_acronym == "CGE" |
                             columns.metadata$structure_acronym == "DTH" |
                             columns.metadata$structure_acronym == "HIP" |
                             columns.metadata$structure_acronym == "LGE" |
                             columns.metadata$structure_acronym == "MGE" |
                             columns.metadata$structure_acronym == "MD" |
                             columns.metadata$structure_acronym == "STR" ] <- "Subcortex"

columns.metadata$Regions[columns.metadata$structure_acronym == "MFC" |
                             columns.metadata$structure_acronym == "DFC" |
                             columns.metadata$structure_acronym == "Ocx" |
                             columns.metadata$structure_acronym == "OFC" |
                             columns.metadata$structure_acronym == "PCx" |
                             columns.metadata$structure_acronym == "TCx" |
                             columns.metadata$structure_acronym == "VFC" |
                             columns.metadata$structure_acronym == "ITC" |
                             columns.metadata$structure_acronym == "STC" |
                             columns.metadata$structure_acronym == "IPC" |
                             columns.metadata$structure_acronym == "V1C"|  
                             columns.metadata$structure_acronym == "M1C" |
                             columns.metadata$structure_acronym == "M1C-S1C" |
                             columns.metadata$structure_acronym == "S1C" |
                             columns.metadata$structure_acronym == "A1C"] <- "Cortex"

columns.metadata$Regions[columns.metadata$structure_acronym == "CBC" |
                             columns.metadata$structure_acronym == "CB" |
                             columns.metadata$structure_acronym == "URL"] <- "Cerebellum"
table(columns.metadata$Regions)

## Conversion to numeric age 
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

## Adding sample name 
columns.metadata <- columns.metadata %>% 
    dplyr::mutate(SampleID = paste(donor_id, age, structure_acronym, stage, sep = "_"))

columns.metadata %<>% 
    as.data.frame() %>%
    dplyr::select("SampleID", everything())

columns.metadata %<>% 
  as.data.frame() %>%
  dplyr::select("SampleID", everything())




## Adding mRIN
columns.metadata$age_for_mRIN <-  gsub("_","", columns.metadata$age)
columns.metadata$donor_name <-  gsub("\\.","_", columns.metadata$donor_name)
columns.metadata = columns.metadata %>% 
    dplyr::mutate(mRIN_id = paste(donor_name, age_for_mRIN,gender, structure_acronym, sep = "//"))
head(columns.metadata$mRIN_id)


bspan.mRIN = readxl::read_xlsx("../BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/mRIN/ncomms8816-s2.xlsx", skip = 1) %>% 
    as.data.frame()
head(bspan.mRIN)

bspan.mRIN %<>% dplyr::select("sample name", "mRIN", "z-score", "P-value")

columns.metadata = merge(columns.metadata, bspan.mRIN, by.x="mRIN_id", by.y = "sample name", all = TRUE, no.dups = TRUE)
columns.metadata = columns.metadata[!is.na(columns.metadata$SampleID),]

columns.metadata = columns.metadata[order(as.numeric(as.character(columns.metadata$column_num))),]
columns.metadata %<>% dplyr::select(-mRIN_id)

write.csv(columns.metadata, "../BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formaatted/BrainSpan-metadata.csv")

## expr 

colnames(counts.matrix) <- columns.metadata$SampleID
rownames(counts.matrix) <- rows.metadata$ensembl_gene_id
counts.matrix %<>% rownames_to_column("EnsemblID")

write.csv(counts.matrix, "../BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formaatted/BrainSpan-exp.csv")
