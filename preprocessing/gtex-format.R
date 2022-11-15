#### 
library(dplyr)

#### Retrieve files from GTEx portal 
#### Save in downloaded data 
#### log and say when the latest download happened 


## Load all gtex files
gtex.raw <- read.delim("gtex_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip=2)

head(gtex.raw)[1:10]
colnames(gtex.raw) <- gsub("\\.", "_", colnames(gtex.raw))
gtex.raw %<>% column_to_rownames("Name")
gtex.brain <- gtex.raw %>% 
  dplyr::select(gtex.md$SampleID)

gtex.raw                           

gtex.md$SampleID

gtex.brain[1:10]
rownames(gtex.brain) <- sub("\\.[0-9]*$", "", rownames(gtex.brain))

write.csv(gtex.brain, "datasets/FormattedData/Gtex/GTEx-exp.csv")
head(gtex.md)






gtex.md.attributes <- read.table("gtex_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
                                 sep = "\t", header=TRUE, fill=TRUE, quote="")
gtex.md.pheno <- read.table("gtex_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", 
                            sep="\t", header=TRUE, fill=TRUE, quote="")

## Format files 
### Fix colnames for the gtex matrix 
colnames(gtex.raw) <- gsub("\\.", "-", colnames(gtex.raw))
gtex.raw %<>% column_to_rownames("Name")

head(gtex.raw)

length(intersect(colnames(gtex.raw), gtex.md.attributes$SAMPID))

#### combine the sample attributes and phenotype attributes files 
x<- transpose(strsplit(gtex.md.attributes$SAMPID, split="-", fixed=TRUE))
SampleName=paste(x[[1]], x[[2]], sep="-")
gtex.md.attributes$SUBJID <- SampleName
pheno <- gtex.md.pheno[match(SampleName, gtex.md.pheno$SUBJID),]
gtex.md.attributes %<>% 
  inner_join(pheno, by = "SUBJID") %>% distinct()



gtex.md.attributes

### Subset brain only datasets from gtex metadata
###### subset samples that belong to the brain and format sample attribute files
brain.samples <- gtex.md.attributes %>% 
  dplyr::filter(grepl("Brain", SMTS)) 
colnames(brain.samples)[1] <- c("SampleID")
head(brain.samples,10)[1:10]  ### 3326 samples 
table(brain.samples$SMTS)

gtex.brain.exp <- gtex.raw %>% 
  dplyr::select(contains(brain.samples$SampleID)) #### 2642 samples subsetted
dim(gtex.brain.exp)

length(intersect(colnames(gtex.brain.exp), brain.samples$SampleID))

table(is.na(match(colnames(gtex.brain.exp), brain.samples$SampleID)))
dim(brain.samples)

length(unique(brain.samples$SampleID))

#### match metadata to exp matrix 

brain.samples <- brain.samples[which(brain.samples$SampleID %in% colnames(gtex.brain.exp)),]

table(is.na(match(colnames(gtex.brain.exp), brain.samples$SampleID)))

rownames(gtex.brain.exp) <- gsub('\\.[0-9]*$', '', rownames(gtex.brain.exp))
gtex.brain.exp %<>% 
  rownames_to_column("EnsemblID")


### changing the names of the columns in gtex
sample.names <- read.xls("gtex_data/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx", header=TRUE)
sample.names <- sample.names[1:2]
sample.names
pheno <- read.xls("gtex_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx", header=TRUE)
pheno <- pheno[1:2]

pheno
colnames(brain.samples)[1] <- c("SAMPID")

gtex.names <- rbind(sample.names, pheno)
head(gtex.names, 20)


brain.samples %<>% 
  dplyr::select(gtex.names$VARNAME)

brain.samples

### load formatted names and save them as a file 
gtex.formatted.names <- read.csv("datasets/FormattedData/Gtex/GTEx-metadata.csv", row.names = 1)
gtex.formatted.names <- colnames(gtex.formatted.names) %>% as.data.frame()
gtex.formatted.names <- gtex.formatted.names[c(-(68:71)),]
gtex.formatted.names

colnames(brain.samples) <- gtex.formatted.names

brain.samples

#### structure acronym 
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Amygdala"] <- "AMY"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Cerebellum"] <- "CB"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Hypothalamus"] <- "HYP"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Substantia nigra"] <- "SNA"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Anterior cingulate cortex (BA24)"] <- "ACC"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Cortex"] <- "CTX"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Nucleus accumbens (basal ganglia)"] <- "NAC"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Caudate (basal ganglia)"] <- "CGE"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Frontal Cortex (BA9)"] <- "DLPFC"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Putamen (basal ganglia)"] <- "PUT"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Cerebellar Hemisphere"] <- "CBC"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Hippocampus"] <- "HIP"
brain.samples$StructureAcronym[brain.samples$Structure == "Brain - Spinal cord (cervical c-1)"] <- "SC1"

#### add regions 
brain.samples$Regions[brain.samples$StructureAcronym == "SC1"] <- "SpinalCord"
brain.samples$Regions[brain.samples$StructureAcronym == "AMY" |
                        brain.samples$StructureAcronym == "CGE" |
                        brain.samples$StructureAcronym == "SNA" |
                        brain.samples$StructureAcronym == "PUT" |
                        brain.samples$StructureAcronym == "NAC" |
                        brain.samples$StructureAcronym == "HYP" |
                        brain.samples$StructureAcronym == "HIP" ] <- "Subcortex"

brain.samples$Regions[brain.samples$StructureAcronym == "ACC" | 
                        brain.samples$StructureAcronym == "DLPFC" |
                        brain.samples$StructureAcronym == "CTX"] <- "Cortex"


brain.samples$Regions[brain.samples$StructureAcronym == "CB" |
                        brain.samples$StructureAcronym == "CBC"] <- "Cerebellum"


##### 
brain.samples$Period <- "Postnatal"
brain.samples$Diagnosis <- "Control"
brain.samples$AgeInterval <- paste(brain.samples$AgeInterval, "yrs", sep = "")
#####

brain.samples$Sex <- gsub("1", "M", brain.samples$Sex)
brain.samples$Sex <- gsub("2", "F", brain.samples$Sex)

brain.samples
### save the metadata and gene expression files 
write.csv(gtex.formatted.names, "colnames_gtex.csv")
write.csv(gtex.brain.exp, "datasets/FormattedData/Gtex/GTEx-exp.csv")
write.csv(brain.samples, "datasets/FormattedData/Gtex/GTEx-metadata.csv")
