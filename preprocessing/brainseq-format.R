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

## Weighting technical replicates
path.data <- file.path("Sanity_check/Raw_data/DownloadedData/BrainSeq/")
## Load
load(file.path(path.data, "rse_gene_unfiltered.Rdata"))
x <- rse_gene@colData
x <- as.data.frame(x)
x <- as.data.frame(t(x))

## Get indices
replicated <- colnames(x)[grep(",", x["SAMPLE_ID",])]

## Subset
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

#write.csv(w, "/Volumes/Data1/PROJECTS/Urwah/integration_datasets/FormattedData/BrainSeq/brainseq_metadata_new_weighted.csv")

## Adding composition data
# load file F9 from https://github.com/LieberInstitute/brainseq_phase2#public-files
load(file.path(path.data, "methprop_pd.Rdata"))

# rename
comp <- as.data.frame(pd)

# filter to composition columns only
comp <- comp[,57:64]

# sanity check
rowSums(comp) # equals 1, per the sum-to-one constrain in deconvolution

# confirm row names are the same
m <- match(rownames(comp), rownames(w)) # they are
final <- cbind(w, comp[m,])

#write.csv(final, "/Volumes/Data1/PROJECTS/Urwah/integration_datasets/FormattedData/BrainSeq/brainseq_metadata_new_weighted_deconv.csv")

### BrainSeq 
### Here we will be formatting the BrainSeq data by the Lieber Institute 
### It is important to note that deconvolution values have been added to the original BrainSeq data in its raw format 
### using the "Weighting Technical Replicates & Adding Composition.R" script. The same script has also been used to fix character variables
### As a result, we are formatting the file that contains the deconvolved values
### of note: the file produced by the previous script is known as brainseq_metadata_new_weighted_deconv.csv

### Load file 
#path.data <- file.path()
#bseq.md <- read.csv(file.path(path.data, "brainseq_metadata_new_weighted_deconv.csv"))



### Adding missing data and quick fixes 

### Period 
final$Period <- ifelse(final$Age > 0, "Postnatal", "Prenatal")

### Age Intervals 
final %<>% mutate(AgeInterval =
                            case_when(
                              between(Age, 1, 1.5) ~ "0-5mos", 
                              between(Age,2,5) ~ "19mos-5yrs", 
                              between(Age,6,11) ~ "6-11yrs", 
                              between(Age, 12,19) ~ "12-19yrs",
                              between(Age, 20,39) ~ " 20-29yrs",
                              between(Age, 30,39) ~ "30-39yrs",
                              between(Age, 40, 49) ~ "40-49yrs", 
                              between(Age, 50, 59) ~ "50-59yrs", 
                              between(Age, 60, 69) ~ "60-69yrs", 
                              between(Age, 70, 79) ~ "70-79yrs", 
                              between(Age, 80, 89) ~ "80-89yrs", 
                              between(Age,90,99) ~ "90-99yrs"))


### Luckily base fuction findInterval works here 
final$AgeInterval[findInterval(final$Age, c(-.7, -.47))] <- "4-7pcw"
final$AgeInterval[findInterval(final$Age, c(-0.46, -0.41))] <- "8-9pcw"
final$AgeInterval[findInterval(final$Age), c(-0.38, -0.33)] <- "10-12pcw"


table(final$Age)
### Regions 
final$Regions[final$Region == "DLPFC"] <- "Cortex"
final$Regions[final$Region == "HIPPO"] <- "Subcortex"
final$Region <- gsub("HIPPO", "HIP", final$Region)


### Name 


bseq.bithub <- read.csv("Sanity_check/BITHub_formatted_data/FormattedData/BrainSeq/BrainSeq-metadata.csv")
table(bseq.bithub$StructureAcronym, bseq.bithub$Regions)

head(bseq.bithub,10)
head(final,10)
