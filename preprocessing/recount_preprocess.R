### Redoing annotations for each of the new datasets 
library(recount3)
library(DT)
library(dplyr)
library(tidyr)
source("functions.R")
library(pander)
library(gridExtra)
library(variancePartition)
library(corrplot)

recount_annot = read.csv("../annotations/Recount-annotations.csv", header=TRUE, check.names = FALSE, sep = ",")
hbdr_annot =read.csv("../annotations/HDBR-metadata.annot.csv", header=TRUE)

hbdr_annot= hbdr_annot[1:15,]

hdbr = recount3::create_rse_manual(
  project = "ERP016243",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

colData(hdbr)

metadata(hdbr)
hbdr_expanded = expand_sra_attributes(hdbr)

### from tutorial 
### Formatting metadata and adding new information in 
#md.hdbr = read.table("/home/neuro/Documents/BrainData/Bulk/HDBR/E-MTAB-4840-experiment-design.tsv", 
 #                    fill=TRUE, sep = "\t", header=TRUE )


#save(md.hdbr, file = "/home/neuro/Documents/Brain_integrative_transcriptome/brain_transcriptome/annotations/HDBR.Rda")

load("../annotations/HDBR.Rda")

md.hdbr = md.hdbr[,1:15]

colnames(md.hdbr) = hbdr_annot$BITColumnName[match(colnames(md.hdbr), hbdr_annot$OriginalMetadataColumnName)]

md.hdbr %<>% 
  as.data.frame() %>% 
  dplyr::select(-c(Block, OntologyIndividual,  KaryotypeOntology, Organism, OrganismOntology, 
                   StructureOntology, HemisphereOntology, OntologyAge)) %>% 
  mutate(StructureAcronym =  add_feature(.$Structure, structure_acronym ), 
         AgeInterval = add_feature(.$Age, age_intervals)) %>% as.data.frame() #%>%
  mutate(Regions = add_feature(.$StructureAcronym, regions))


add_feature = function(feature_column, features){
    as.vector(sapply(feature_column, function(x){
      names(features)[sapply(features, function(f) x %in% f)]})) 
  }
  
  
class(y$StructureAcronym) 
md.hdbr= data.frame(md.hdbr)

x <- as.list(md.hdbr)  
(y <- do.call(cbind, x))
 
y= as.data.frame(y)  
md.hdbr$Sample.Characteristic.developmental.stage.


metadata = colData(hdbr)

table(md.hdbr$Sample.Characteristic.organism.part.)
metadata %>% 
  as.data.frame() %>% 
  dplyr::select(contains("sra")) %>% 
  DT::datatable()


####


### Ramakar et al (2017)

#### Gene level data 
ramakar = recount3::create_rse_manual(
  project = "SRP073813",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

### Format SRA data
sra_data = colData(ramakar) %>% 
  as.data.frame() %>% 
  dplyr::select(contains("sra")) 


md.sra = sra_data %>%
  dplyr::select("sra.sample_attributes") %>%
  separate(sra.sample_attributes, c("Age", "pH", "StructureAcronym", "Diagnosis", "Ethnicity", "Sex",
                                    "PMI", "SourceName", "Tissue"), sep="\\|") %>% 
  mutate(Age = gsub(".*;;","", Age), 
         pH = gsub(".*;;","", pH), 
         StructureAcronym = gsub(".*;;","", StructureAcronym),
         Diagnosis = gsub(".*;;","", Diagnosis), 
         Ethnicity = gsub(".*;;","", Ethnicity), 
         Sex = gsub(".*;;","", Sex), 
         PMI = gsub(".*;;","", PMI)) %>% 
  dplyr::select(-c("SourceName", "Tissue")) %>% 
  mutate(Sex = ifelse(Sex == "male", "M", "F"), 
         StructureAcronym = gsub("nAcc","NAC", StructureAcronym), 
         StructureAcronym = gsub("AnCg", "ACC", StructureAcronym)) %>%
  mutate(Regions =  add_feature(.$StructureAcronym, regions),
        Age_rounded = as.character(sapply(na.omit(as.numeric(.$Age)), num_to_round)), 
        Period = c("Postnatal")) %>%
  mutate(AgeInterval = add_feature(.$Age_rounded , age_intervals)) %>%
  rownames_to_column("SampleID")

### format recount annotations 

md = colData(ramakar) %>% 
  as.data.frame() %>%
  dplyr::select(-contains("2"))
colnames(md) = recount_annot$BITHubName[match(colnames(md), recount_annot$RecountName)]
md = full_join(md, md.sra)


## Removing columns with NA values
md.clean = md %>% select_if(~ !any(is.na(.)))
md.clean = md.clean[vapply(md.clean, function(x) length(unique(x)) > 1, logical(1L))]
### Check correlations of samples 


M = cor(data.matrix(md.clean), use = "complete.obs")
corrplot(M, method = 'number', order='AOE')


### FANTOM5

f5 =read.table("/home/neuro/Documents/BrainData/Bulk/Fantom5/hg19.cage_peak_phase1and2combined_ann.txt", fill = TRUE)
head(f5)
