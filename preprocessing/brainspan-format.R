suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tools)
  library(magrittr)
  library(tibble)
  library(plyr)
  library(pheatmap)
  library(readxl)
  library(gdata)
  library(biomaRt)
  library(data.table)
  library(pander)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  library(stringr)
  library(SummarizedExperiment)
  library(kableExtra)
  library(scales)
  library(viridis)})

## Remove existing lists 
rm(list=ls())

## BrainSpan 

### Download new data from BrainSpan 
#download.file('https://www.brainspan.org/api/v2/well_known_file_download/267666525', 
#destfile = "DownloadedData/BrainSpan", method = "wget", extra = "-r -p --random-wait")


### DEFINE PATHS HERE 
path.data <- file.path("Sanity_check/Raw_data/DownloadedData/genes_matrix_csv/") ## Should be changed to RNA 
out.path <- file.path("")

## Load all BrainSpan files
columns.metadata <- read.csv(file.path(path.data, "columns_metadata.csv"), header = TRUE)
counts.matrix <- read.csv(file.path(path.data, "expression_matrix.csv"), header= FALSE, row.names= 1)
rows.metadata <- read.csv(file.path(path.data, "rows_metadata.csv"))


### Fix format and levels for age
columns.metadata$age <-  gsub(" ","-", columns.metadata$age)
levels(columns.metadata$age) <- c("8-pcw","12-pcw", "13-pcw", "16-pcw", 
                                  "17-pcw", "19-pcw", "21-pcw",
                                  "24-pcw", "25-pcw", "26-pcw", "35-pcw", 
                                  "37-pcw", "4-mos","6-mos", "10-mos", "1-yrs", 
                                  "2-yrs","3-yrs", "4-yrs", "8-yrs", "11-yrs", 
                                  "13-yrs", "15-yrs", "18-yrs", "19-yrs", "21-yrs", 
                                  "23-yrs", "30-yrs", "36-yrs", "37-yrs", "40-yrs")



## ADDING STAGES 
## Format all missing information in the metadata file using the Technical White paper 
## On top of adding stages, we have further subsetted stages for adult stages 
## e.g according to the technical white paper, the adult ages of
## 20yrs-60yrs would be stage 11 (s11), however we have binned them by decades  

columns.metadata$stage[columns.metadata$age == "8-pcw" | columns.metadata$age == "9-pcw"] <- "s02a" 
columns.metadata$stage[columns.metadata$age == "12-pcw"] <- "s02b"
columns.metadata$stage[columns.metadata$age == "13-pcw"] <- "s03a"
columns.metadata$stage[columns.metadata$age == "16-pcw" | columns.metadata$age == "17-pcw"] <-"s03b"
columns.metadata$stage[columns.metadata$age == "19-pcw" | columns.metadata$age == "21-pcw" | columns.metadata$age == "24-pcw"] <- "s04"
columns.metadata$stage[columns.metadata$age == "25-pcw" | columns.metadata$age == "26-pcw" | columns.metadata$age == "35-pcw" | columns.metadata$age == "37-pcw"] <- "s05"
columns.metadata$stage[columns.metadata$age == "4-mos" | columns.metadata$age == "6-mos"] <- "s06"
columns.metadata$stage[columns.metadata$age == "10-mos" | columns.metadata$age == "1-yrs"] <- "s07"
columns.metadata$stage[columns.metadata$age == "2-yrs" | columns.metadata$age == "3-yrs" | columns.metadata$age == "4-yrs"] <- "s08"
columns.metadata$stage[columns.metadata$age == "8-yrs"| columns.metadata$age == "11-yrs"] <- "s09"
columns.metadata$stage[columns.metadata$age == "13-yrs" | columns.metadata$age == "15-yrs" | columns.metadata$age == "18-yrs" | columns.metadata$age == "19-yrs"] <-  "s10"
columns.metadata$stage[columns.metadata$age == "21-yrs" |
                         columns.metadata$age == "23-yrs"] <- "s11"
columns.metadata$stage[columns.metadata$age == "30-yrs" |
                         columns.metadata$age == "36-yrs" | columns.metadata$age == "37-yrs"] <- "s12"
columns.metadata$stage[columns.metadata$age == "40-yrs"] <- "s13"


levels(columns.metadata$stage) <- c("s02a", "s02b", "s03a", "s03b", "s04", "s05", "s06", 
                                    "s07", "s08", "s09", "s10", "s11", "s12", "s13")


## ADDING PERIOD 

columns.metadata$Period[columns.metadata$stage == "s02a"|
                          columns.metadata$stage == "s02b"|
                          columns.metadata$stage == "s03a"|
                          columns.metadata$stage == "s03b"|
                          columns.metadata$stage == "s04"|
                          columns.metadata$stage == "s05"] <- "Prenatal"

columns.metadata$Period[columns.metadata$stage == "s06"|
                          columns.metadata$stage == "s07"|
                          columns.metadata$stage == "s08"|
                          columns.metadata$stage == "s09"|
                          columns.metadata$stage == "s10"|
                          columns.metadata$stage == "s11" |
                          columns.metadata$stage == "s12" |
                          columns.metadata$stage == "s13"] <- "Postnatal"




### Addition of Regions 
#### Regions indicate the wider area of the brain a structure 
#### is a part of 

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



## Conversion to numeric values for ages 


### conversion of pcw to nuermic ages -- using the formula -((40-x) / 52)
columns.metadata$AgeNumeric[grepl("pcw", columns.metadata$age, ignore.case = TRUE)]<-
  columns.metadata$age[grepl("pcw", columns.metadata$age)] %>%  str_remove("-pcw")%>% 
  as.numeric() %>% `-` (40) %>% divide_by(52) ### this is the correct way of doing things 

#bspan.metadata$AgeNumeric[grepl("pcw", bspan.metadata$Age, ignore.case = TRUE)] <- 
 # bspan.metadata$Age[grepl("pcw",  bspan.metadata$Age)] %>%  str_remove("-pcw") %>% 
  #as.numeric() %>% -(40 - .) %>% divide_by(52) ### Method used originally - does not work 


#### conversion of mos to numeric ages 
columns.metadata$AgeNumeric[grepl("mos", columns.metadata$age, ignore.case = TRUE)] <- 
  columns.metadata$age[grepl("-mos", columns.metadata$age)] %>%  str_remove("-mos") %>% 
  as.numeric() %>% divide_by(12)


#### conversion of yrs to numeric ages - already numeric, just need to remove the -yrs suffix
columns.metadata$AgeNumeric[grepl("yrs", columns.metadata$age, ignore.case = TRUE)] <- 
  columns.metadata$age[grepl("-yrs", columns.metadata$age)] %>%  str_remove("-yrs") %>% 
  as.numeric 



## Add diagnosis 
### According to the technical white paper, all BrainSpan data are controls 

columns.metadata$Diagnosis <- "Control"


## Add age interval 
columns.metadata$AgeInterval[grep("s02a", columns.metadata$stage)] <- "8-9pcw"
columns.metadata$AgeInterval[grepl("s02b", columns.metadata$stage, ignore.case = TRUE)] <- "10-12pcw"
columns.metadata$AgeInterval[grepl("s03a", columns.metadata$stage, ignore.case = TRUE)] <- "13-15pcw"
columns.metadata$AgeInterval[grepl("s03b", columns.metadata$stage, ignore.case = TRUE)] <- "16-18pcw"
columns.metadata$AgeInterval[grepl("s04", columns.metadata$stage, ignore.case = TRUE)] <- "19-24pcw"
columns.metadata$AgeInterval[grepl("s05", columns.metadata$stage, ignore.case = TRUE)] <- "25-38pcw"
columns.metadata$AgeInterval[grepl("s06", columns.metadata$stage, ignore.case = TRUE)] <- "0-5mos"
columns.metadata$AgeInterval[grepl("s07", columns.metadata$stage, ignore.case = TRUE)] <- "6-18mos"
columns.metadata$AgeInterval[grepl("s08", columns.metadata$stage, ignore.case = TRUE)] <- "19mos-5yrs"
columns.metadata$AgeInterval[grepl("s09", columns.metadata$stage, ignore.case = TRUE)] <- "6-11yrs"
columns.metadata$AgeInterval[grepl("s10", columns.metadata$stage, ignore.case = TRUE)] <- "12-19yrs"
columns.metadata$AgeInterval[grepl("s11", columns.metadata$stage, ignore.case = TRUE)] <- "20-29yrs"
columns.metadata$AgeInterval[grepl("s12", columns.metadata$stage, ignore.case = TRUE)] <- "30-39yrs"
columns.metadata$AgeInterval[grepl("s13", columns.metadata$stage, ignore.case = TRUE)] <- "40-49yrs"





### Add sampleID names 
columns.metadata <- columns.metadata %>% 
  dplyr::mutate(SampleID = paste(donor_id, age, structure_acronym, stage, sep = "_"))

columns.metadata %<>% 
  as.data.frame() %>%
  dplyr::select("SampleID", everything())

table(columns.metadata$AgeInterval, columns.metadata$AgeNumeric)

###Change column names - we edited colnames for BIThub for
### consistency across datasets in a seperate file
### now we will using these to change the colnames in the dataset before saving  



bspan.names <- read.csv("Sanity_check/Raw_data/DownloadedData/BrainSpan-metadata-annot.csv")
bspan.names
colnames(columns.metadata)[c(1:10)] <-bspan.names$BITColumnName[c(1:10)]
colnames(columns.metadata)[8] <- "StructureAcronym"

columns.metadata
### Add mRIN values to BrainSpan 
bs.mRIN <- read.csv("../CP_gene_review/BrainSpan-mRIN.csv")[,-1]
columns.metadata = cbind(columns.metadata, bs.mRIN)
columns.metadata

columns.metadata %<>% 
  as.data.frame() %<>%
  dplyr::select(-c("name_for_matrix", "mRIN_name", "sample.name","z.score", "P.value"))

###  Add SampleID as colname and geneID as rowname in the brainspan exp file 
colnames(counts.matrix) <- columns.metadata$SampleID
rownames(counts.matrix) <- rows.metadata$ensembl_gene_id
counts.matrix %<>% rownames_to_column("EnsemblID")

### save data

head(bspan.names)
write.csv(columns.metadata, "BrainSpan-metadata.csv")
write.csv(counts.matrix, "BrainSpan-exp.csv")


### Comparison with the current BITHub file 

#bspan.bithub <- read.csv("Sanity_check/BITHub_formatted_data/FormattedData/BrainSpan/BrainSpan-metadata.csv")
#head(bspan.bithub)

#table(is.na(match(bspan.bithub$SampleID, columns.metadata$SampleID))) ### sampleIDs match order 

#table(columns.metadata$Period)
#columns.diff <- subset(columns.metadata, !(AgeNumeric %in% bspan.bithub$AgeNumeric))

#bspan.diff <-  subset(bspan.bithub, !(AgeNumeric %in% columns.metadata$AgeNumeric))
#head(bspan.diff$AgeNumeric,10)
#head(columns.diff$AgeNumeric,10)


### change 
### Previous method 
#columns.metadata$AgeNumeric[grepl("pcw", columns.metadata$Age, ignore.case = TRUE)] <- 
 # columns.metadata$Age[grepl("pcw",  columns.metadata$Age)] %>%  str_remove("-pcw") %>% 
#  as.numeric() %>% -(40 - .) %>% divide_by(52)

### Gavins method 
# step 2: convert PCW into a negative age
#which.pcw <- grep("-pcw", levels(columns.metadata$Age))
#pcw <- levels(columns.metadata$Age)[which.pcw]

#pcw <- gsub("-pcw", "", pcw)
#pcw <- as.numeric(pcw)
#pcw <- -(40 - pcw) # how many weeks before birth
#pcw <- pcw / 52 # converts weeks into a fraction of the year
#unique(pcw)
#unique(columns.metadata$AgeNumeric)
#unique(bspan.bithub$AgeNumeric)
#levels(columns.metadata$AgeNumeric)[which.pcw] <- pcw
#columns.metadata$AgeNumeric
