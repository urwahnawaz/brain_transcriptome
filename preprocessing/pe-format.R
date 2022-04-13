##### Load original metadata
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
  library(splitstackshape)
  library(viridis)})


rm(list=ls())

path.data <- file.path("Sanity_check/Raw_data/DownloadedData/Psychencode")
old <- read.csv(file.path(path.data, "Job-154061626746631166818172835.csv"))

table(old$ageDeath)
old
##### Gavin's script 
## Aside: first, convert age to numeric
# step 1: convert "90+" to "91"
levels(old$ageDeath)[grep("\\+", levels(old$ageDeath))] <- "91"

# step 2: convert PCW into a negative age
which.pcw <- grep("PCW", levels(old$ageDeath))
pcw <- levels(old$ageDeath)[which.pcw]

pcw <- gsub(" PCW", "", pcw)
pcw <- as.numeric(pcw)
pcw <- -(40 - pcw) # how many weeks before birth
pcw <- pcw / 52 # converts weeks into a fraction of the year
levels(old$ageDeath)[which.pcw] <- pcw

# step 3: convert to numeric (also turns a blank "" into NA)
old$ageDeath <- as.numeric(as.character(old$ageDeath))

## Load composition
comp <- read_excel(file.path(path.data,"DER-24_Cell_fractions_Normalized.xlsx"))
comp <- as.data.frame(comp)
rownames(comp) <- comp[,1]
comp <- comp[,-1]
comp <- as.data.frame(t(comp))

## Augment
m <- match(old$individualID, rownames(comp))
new <- cbind(old, comp[m,])

## Save
#write.csv(new, "/Volumes/Data1/PROJECTS/Urwah/integration_datasets/FormattedData/PsychEncode/Job-154061626746631166818172835-UpdatedPerScript.csv")


#### PsychEncode data formatting - Urwah's code
### Here we will be formatting the psychencode data 
### It is important to note that deconvolution values have been added to the original psychencode data in its raw format 
### using the AddingComposition.R script. The same script has also been used to convert ages into numeric ages 
### As a result, we are formatting the file that contains the deconvolved values and the numeric ages 
### of note: the file produced by the previous script is known as Job-154061626746631166818172835-UpdatedPerScript.csv 

#path.data <- file.path("../../CoC/Data/Bulk_data/DownloadedData/genes_matrix_csv/") ### Change directory to rna 


#PE.metadata <- read.csv(file.path(path.data,"deconv_data/Job-154061626746631166818172835-UpdatedPerScript.csv"))
#PE.exp <- read.table(file.path(path.data, "DER-01_PEC_Gene_expression_matrix_normalized.txt"), header=TRUE, row.names = 1)
PE.exp <- read.table("DER-02_PEC_Gene_expression_matrix_TPM.txt", header=TRUE, row.names = 1)
### what is the difference between the two files?? 

### change to DER-02 and ask during lab meeting 

#### We are removing samples that contain BP in their diagnosis 
no.BP <- new %>%
  dplyr::filter(diagnosis == "Affective Disorder" | diagnosis == "Autism Spectrum Disorder" | diagnosis == "Bipolar Disorder" |
                  diagnosis == "Control" | diagnosis == "Schizophrenia")

table(no.BP$diagnosis)



### Adding structure and structure acronym information 
no.BP$Structure <- c("Prefrontal Cortex")
no.BP$StructureAcronym<- c("PFC")

### Pre or post-natal samples 

no.BP$Period <- ifelse(no.BP$ageDeath > 0, "Postnatal", "Prenatal")

### Age intervals 
no.BP %<>% mutate(AgeInterval =
                    case_when(#between(ageDeath,-0.7, -0.47), ~ "4-7pcw",  
                      #between(ageDeath, -0.46, -0.41), ~ "8-9pcw",  ### UNFORTUNATELY case_when does not work
                      #between(ageDeath, -0.39, -0.37), ~ "10-12pcw", ### with negative values 
                      between(ageDeath, 0, 1.5) ~ "0-5mos", 
                      between(ageDeath,2,5.99) ~ "19mos-5yrs", 
                      between(ageDeath,6,11.99) ~ "6-11yrs", 
                      between(ageDeath, 12,19.99) ~ "12-19yrs",
                      between(ageDeath, 20,29.99) ~ " 20-29yrs",
                      between(ageDeath, 30,39.99) ~ "30-39yrs",
                      between(ageDeath, 40, 49.99) ~ "40-49yrs", 
                      between(ageDeath, 50, 59.99) ~ "50-59yrs", 
                      between(ageDeath, 60, 69.99) ~ "60-69yrs", 
                      between(ageDeath, 70, 79.99) ~ "70-79yrs", 
                      between(ageDeath, 80, 89.99) ~ "80-89yrs", 
                      between(ageDeath,90,99.99) ~ "90-99yrs"))


no.BP

### Luckily base fuctions would be of better use here 

no.BP$AgeInterval[no.BP$ageDeath == -0.5945210] <- "8-9pcw"
no.BP$AgeInterval[no.BP$ageDeath >= -0.52 & no.BP$ageDeath <= -0.47] <- "13-15pcw"
no.BP$AgeInterval[no.BP$ageDeath >= -0.47 & no.BP$ageDeath <= -0.42] <- "16-18pcw"
no.BP$AgeInterval[no.BP$ageDeath >= -0.41 & no.BP$ageDeath <= -0.33] <- "19-24pcw"
no.BP$AgeInterval[no.BP$ageDeath >= -0.27 & no.BP$ageDeath <= -0.090] <- "25-38pcw"


no.BP$AgeInterval[findInterval(no.BP$ageOnset, c(-0.52, -0.48))] <- "13-15pcw"
no.BP$AgeInterval[findInterval(no.BP$ageOnset, c(-0.47, -0.42))] <- "16-18pcw"
no.BP$AgeInterval[findInterval(no.BP$ageOnset, c(-0.41, -0.33))] <- "19-24pcw"
no.BP$AgeInterval[findInterval(no.BP$ageOnset, c(-0.2684930, -0.0958900))] <- "25-38pcw"
no.BP$AgeInterval[no.BP$ageDeath == -0.0191780] <- "39-40pcw"

table(no.BP$AgeInterval)
table(no.BP$ageDeath, no.BP$AgeInterval)
### subset samples from the expression matrix
pe.data <- PE.exp %>% 
  dplyr::select(contains(no.BP$individualID)) 

BP.final <- no.BP %>% dplyr::filter(individualID %in% colnames(pe.data))
BP.final <- BP.final[,-1]
BP.final %<>%
  dplyr::select(individualID, everything())

pe.data %<>% rownames_to_column("EnsemblID")




pe.metadata.annotations<- colnames(BP.final) %>% as.data.frame()
length(colnames(BP.final))
colnames(BP.final)
head(BP.final)
colnames(BP.final)[1:19] <- c("SampleID",
                                     "RowID", 
                                     "RowVersion", 
                                     "ContributingStudy", 
                                     "SampleIDSource", 
                                     "Diagnosis", 
                                     "Sex", 
                                     "Ethnicity", 
                                     "AgeNumeric", 
                                     "ageOnset", 
                                     "causeDeath", 
                                     "brainWeight", 
                                     "height", 
                                     "weight", 
                                     "ageBiopsy", 
                                     "smellTestScore",
                                     "smoker", 
                                     "notes",
                                     "Capstone_4")


### compare with the formatted data 
write.csv(BP.final, "PsychEncode-metadata.csv")
write.csv(pe.data, "PsychEncode-exp.csv")

pe.bithub <- read.csv("Sanity_check/BITHub_formatted_data/FormattedData/Psychencode/PsychEncode-metadata.csv")
table(is.na(match(BP.final$SampleID, pe.bithub$SampleID))) ### sample names match
