# Loading all libraries
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


## Weighting technical replicates
path.data <- file.path("../../../BrainData/Bulk/Brainseq/")
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

w

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



### Adding missing data and quick fixes 

### Period 
final$Period <- ifelse(final$Age > 0, "Postnatal", "Prenatal")

### Age Intervals 
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

final$AgeInterval[final$Age == -0.5945210] <- "8-9pcw"
final$AgeInterval[final$Age >= -0.52 & final$Age <= -0.47] <- "13-15pcw"
final$AgeInterval[final$Age >= -0.47 & final$Age <= -0.42] <- "16-18pcw"
final$AgeInterval[final$Age >= -0.41 & final$Age <= -0.33] <- "19-24pcw"
final$AgeInterval[final$Age >= -0.27 & final$Age <= -0.090] <- "25-38pcw"



table(final$Age, final$AgeInterval)
### Regions 
final$Regions[final$Region == "DLPFC"] <- "Cortex"
final$Regions[final$Region == "HIPPO"] <- "Subcortex"
final$Region <- gsub("HIPPO", "HIP", final$Region)

final

### Save the file 

write.csv("../../../BrainData/Bulk/Brainseq/Formatted/", "BrainSeq-metadata.csv")


### changes to on doing age intervals 

### load brainseq
dir = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq")
load(file.path(dir, "rse_gene_unfiltered.Rdata"), envir = .GlobalEnv)



num_to_pcw = function(age){
  round(age * 52 + 40)
}

num_to_mos = function(age){
  round(age * 12 )
}

num_to_mos( 1.623589 )

?round
num_to_pcw(-0.402739726)

w %>% 
  dplyr::filter(Age < 0) %>% 
  mutate(ageNN = paste0(num_to_pcw(Age)))

paste0(c(30:39), " yrs")

table(w$age_pcw)
w$age_pcw[w$Age < 0] = paste0(num_to_pcw(w$Age[w$Age < 0]), " pcw") %>% as.data.frame()
trial =w %>% mutate(age_interval = as.character(cut(Age, seq(-1, 100, by = 10)))) %>%
  mutate(AgeInterval = sapply(age_interval, function(i) {
    paste0( as.numeric(gsub("^\\(([-0-9]+),.+", "\\1", i)) + 1,
            "-", as.numeric(gsub(".+,([0-9]+)\\]$", "\\1", i)), "yrs")})) %>% 
  mutate(Age_Interval = add_feature(.$age_pcw, age_intervals))

table(trial$AgeInterval)

md = md %>% mutate(Period = ifelse(.$Age > 0, "Postnatal", "Prenatal"), 
                   Region = gsub("HIPPO", "HIP", .$Region)) %>%
  mutate(Regions = add_feature(.$Region, regions)) %>% 
  mutate(age_interval = as.character(cut(Age, seq(-1, 100, by = 10)))) %>%
  mutate(AgeInterval = sapply(age_interval, function(i) {
    paste0( as.numeric(gsub("^\\(([-0-9]+),.+", "\\1", i)) + 1,
            "-", as.numeric(gsub(".+,([0-9]+)\\]$", "\\1", i)), "yrs")})) %>% 
  dplyr::select(-age_interval)

mutate(Regions = add_feature(.$StructureAcronym, regions), 




?grepl

outdir = file.path("/home/neuro/Documents/Brain_integrative_transcriptome/Results/Formatted")