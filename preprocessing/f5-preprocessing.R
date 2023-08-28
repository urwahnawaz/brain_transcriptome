### Fantom5 data 
source("functions.R")
library(recount3)
library(magrittr)
library(tibble)
library(reshape2)
library(SummarizedExperiment)
library(corrplot)
library(dplyr)
library(ggvenn)
# 
# f5_data = read.table("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5PHASE1and2FREEZEDPIclustershuman(robustsetwithexpression).osc", header=TRUE)
# f5.md =read.delim("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5PHASE1and2FREEZEDPIclustershuman(robustsetwithexpression)_wMD.csv",
#                   skip =9)
# 
# f5.md = f5.md[,-1]
# f5.md = f5.md[,-3]
# head(f5.md )
# 
# 
# trial = f5.md %>% 
#   mutate(MD_name = sapply(str_split(Metadata, pattern = "\\+AD0-"), function(x) paste(x[1:2], collapse =""))) %>% 
#   mutate(variable = sapply(str_split(MD_name, pattern = "-"), function(x) paste(x[1:2], collapse = ""))) 


### Fantom5 data processing 

## Expression data 
f5_osc =  read.delim("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5/FANTOM5PHASE1and2gencodev19transcriptexpression(robustpromoter-50bp5primeend).osc")
names_in_osc = colnames(f5_osc)
names_in_osc %<>% as.data.frame() %>% set_colnames("Colnames") %>% 
  mutate(SampleID = gsub("exp.rle.", "", Colnames)) %>% 
  mutate(SampleID = gsub("_ctss", "", SampleID))
names_in_osc[-(1:7),]
names_in_osc = names_in_osc[-(1:7),]

### Data from Hamid 
f5_hamid = read.delim("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5/Fantom5SamplesAnnotated.csv", sep = ",")

f5_hamid %<>% 
  mutate(SampleID = gsub("\\..*", "", Library))


x = list( "Hamid's MD" = f5_hamid$SampleID, 
          "F5 OSC colnames" = names_in_osc$SampleID)

venn= ggvenn(x[c(1,2)], fill_color = c("#2F124B", "#EA2A5F"),
                   set_name_color = "black", text_color = "white", text_size = 6, 
                   stroke_size = 0.5, fill_alpha = 0.6) 


#### Downloaded from CAGE paper 
library(readxl)

HS = read_excel("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5/HumanSamples2.0.sdrf.xlsx") %>% as.data.frame()
colnames(HS)[2] = c("FFID")

other_supp = read.table("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5/DRA002747.txt", fill = TRUE, header=TRUE)
head(other_supp)


merged = merge(HS , other_supp, by= "FFID")
colnames(merged) = c("FFID", 
                     "SourceName",
                     "Characteristic", 
                     "Catalog_ID",  
                     "SampleType",
                     "Species",
                     "Sex", 
                     "AgeNumeric",
                     "Stage", 
                     "Tissue", 
                     "CellLot", 
                     "CellType", 
                     "Catalogue ID",
                     "Collaboration", 
                     "Provider", 
                     "Protocol",
                     "ExtrctName",
                     "MaterialType",
                     "LibraryID", 
                     "DRA_sample_accession", 
                     "DRA_experiment_accession",
                     "DRA_run_accession", 
                     "DRA_bam_accession",
                     "DRA_ctss_accession",
                     "Experimenet_method")


x = list( "HS" = HS$FFID, 
          "Other supp" = other_supp$FFID)

ggvenn(x[c(1,2)], fill_color = c("#2F124B", "#EA2A5F"),
       set_name_color = "black", text_color = "white", text_size = 6, 
       stroke_size = 0.5, fill_alpha = 0.6) 


### Another one 

more_supp =  read.table("/home/neuro/Documents/BrainData/Bulk/Fantom5/FANTOM5/41597_2017_BFsdata2017112_MOESM152_ESM/a_assay_Kawaji.txt", fill = TRUE, header=TRUE)


x = list( "Hamid's MD" = f5_hamid$SampleID, 
          "Supp" = other_supp$LibraryID)
ggvenn(x[c(1,2)], fill_color = c("#2F124B", "#EA2A5F"),
       set_name_color = "black", text_color = "white", text_size = 6, 
       stroke_size = 0.5, fill_alpha = 0.6)  ## Doesn't merge with Hamids 


## with colnames 
x = list( "F5 OSC" = names_in_osc$SampleID, 
          "Supp" = merged$LibraryID)
ggvenn(x[c(1,2)], fill_color = c("#2F124B", "#EA2A5F"),
       set_name_color = "black", text_color = "white", text_size = 6, 
       stroke_size = 0.5, fill_alpha = 0.6)  



### other Supp and Hamid's md 
merged %>% head(5)
merged_formatted = merged %>% 
  dplyr::select(-c("Provider", "Catalogue ID", "Stage", 
                   "CellLot", "Catalog_ID")) %>%
  mutate(Sex = gsub("female", "F", Sex),
         Sex = gsub("male", "M", Sex), 
         Sex = gsub("n/a", "NA", Sex),
         Sex = gsub("UNDEFINED_SEX_TYPE", "Undefined or Unknown", Sex), 
         Sex =gsub("unknown", "Undefined or Unknown", Sex),
         Tissue = str_to_sentence(Tissue),
         CellType =str_to_sentence(CellType), 
         DonorID = ifelse(grepl("donor",  Characteristic), Characteristic, NA)) %>% 
  mutate(DonorID = gsub(".*, ", "", DonorID)) %>% 
  mutate(DonorID = str_to_sentence(DonorID),
         CellLine = ifelse(grepl("cell line", Characteristic), Characteristic, NA)) %>% 
  dplyr::rename("SampleID" = "LibraryID")


f5_hamid
