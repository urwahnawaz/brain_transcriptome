## What columns are in the processed GTEx file 
library(dplyr)
library(purrr)
gtex.processed = read.csv("../../../BrainData/Bulk/FormattedData/Gtex/GTEx-metadata.csv", header= TRUE, row.names = 1)



gtex.md.attributes = read.table("../../../BrainData/Bulk/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
                                sep = "\t", header=TRUE, fill=TRUE, quote="")
gtex.md.pheno = read.table("../../../BrainData/Bulk/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", 
                           sep="\t", header=TRUE, fill=TRUE, quote="")

gtex.annot = read.csv("../../../BrainData/Bulk/FormattedData/Gtex/GTEx-metadata-annot.csv", header = TRUE)



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

## Get only brain samples 
gtex.md.attributes %<>% 
  dplyr::filter(grepl("Brain", Structure)) 

## Add other annotations 
gtex.md.attributes$StructureAcronym = as.factor(sapply(gtex.md.attributes$Structure, function(x){
  names(structure_acronym)[sapply(structure_acronym, function(structure) x %in% structure)]})) %>% as.data.frame()

dim(gtex.md.attributes$StructureAcronym)
table(gtex.md.attributes$StructureAcronym)
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Amygdala"] <- "AMY"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Cerebellum"] <- "CB"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Hypothalamus"] <- "HYP"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Substantia nigra"] <- "SNA"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Anterior cingulate cortex (BA24)"] <- "ACC"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Cortex"] <- "CTX"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Nucleus accumbens (basal ganglia)"] <- "NAC"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Caudate (basal ganglia)"] <- "CGE"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Frontal Cortex (BA9)"] <- "DLPFC"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Putamen (basal ganglia)"] <- "PUT"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Cerebellar Hemisphere"] <- "CBC"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Hippocampus"] <- "HIP"
gtex.md.attributes$StructureAcronym[gtex.md.attributes$Structure == "Brain - Spinal cord (cervical c-1)"] <- "SC1"



gtex.md.attributes$Regions = as.factor(sapply(gtex.md.attributes$StructureAcronym, function(x){
  names(regions)[sapply(regions, function(region) x %in% region)]}))


gtex.md.attributes$Period <- "Postnatal"
gtex.md.attributes$Diagnosis <- "Control"
gtex.md.attributes$AgeInterval <- paste(gtex.md.attributes$AgeInterval, "yrs", sep = "")


gtex.md.attributes$Sex <- gsub("1", "M", gtex.md.attributes$Sex)
gtex.md.attributes$Sex <- gsub("2", "F", gtex.md.attributes$Sex)

head(gtex.md.attributes)
