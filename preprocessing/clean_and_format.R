source("functions.R")


# Preprocessing of BrainSpan files 

## In the bspan directory, you must ensure you have the following files 
## columns_metadata.csv
## rows_metadata.csv
## expression_matrix.csv

bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/")

clean_and_format(bspandir,"BrainSpan", outdir) # This function will save the files in the outdir



## Brainspan check code 
# Loading all files 

dir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")
annot = read.csv("../annotations/BrainSpan-metadata-annot.csv")

columns.metadata = read.csv(file.path(dir, "columns_metadata.csv"), header = TRUE)
exp = read.csv(file.path(dir, "expression_matrix.csv"), header= FALSE, row.names= 1)
rows.metadata = read.csv(file.path(dir, "rows_metadata.csv"))

bspan.mRIN = readxl::read_xlsx("../annotations/ncomms8816-s2.xlsx", skip = 1) %>% 
  as.data.frame() %>% dplyr::select("sample name", "mRIN", "z-score", "P-value")

message("Changing BrainSpan column names to BITHub column names")
colnames(columns.metadata) = annot$BITColumnName[match(colnames(columns.metadata), annot$OriginalMetadataColumnName)]

message("Adding additional metadata information")

md = columns.metadata %>% 
  mutate(Stage = add_feature(.$Age, stages), 
         Regions = add_feature(.$StructureAcronym, regions), 
         AgeInterval = add_feature(.$Age, age_intervals), 
         Diagnosis = "Control", 
         Age = gsub(" ","_", .$Age)) %>%  
  mutate(SampleID = paste(DonorID, Age, StructureAcronym, Stage, sep = "_"), 
         age_for_mRIN = gsub("_", "", .$Age), 
         DonorName = gsub("\\.","_", .$DonorName)) %>% 
  mutate('sample name' = paste(DonorName, age_for_mRIN,Sex ,StructureAcronym, sep = "//")) %>% 
  left_join(bspan.mRIN) %>% as.data.frame() %>% 
  dplyr::select("SampleID", everything())

## Adding Age Numeric 
md$AgeNumeric[grepl("pcw", md$Age, ignore.case = TRUE)]<-
  md$Age[grepl("pcw", md$Age)] %>%  str_remove("_pcw")%>% 
  as.numeric() %>% `-` (40) %>% divide_by(52)
md$AgeNumeric[grepl("mos", md$Age, ignore.case = TRUE)] <- 
  md$Age[grepl("_mos", md$Age)] %>%  str_remove("_mos") %>% 
  as.numeric() %>% divide_by(12)
md$AgeNumeric[grepl("yrs", md$Age, ignore.case = TRUE)] <- 
  md$Age[grepl("_yrs", md$Age)] %>%  str_remove("_yrs") %>% 
  as.numeric 

md %<>% mutate(Period = ifelse(.$AgeNumeric >= 0, "Postnatal", "Prenatal"), 
               colname = paste(DonorName, Age, StructureAcronym, sep = "_"))


## Add from additional excel file retrived from Allen Brain Atlas  
md.excel = read_excel(file.path(dir, "BrainSpan-additional.xlsx"),sheet =2, col_names = TRUE, skip =1) %>% 
  as.data.frame() %>% mutate_at(.vars = "AllenInstituteID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>%
  mutate_at(.vars = "Age", .funs = gsub, pattern = "PCW", replacement = "_pcw") %>% 
  mutate_at(.vars = "Age", .funs = gsub, pattern = "M", replacement = "_mos") %>%
  mutate_at(.vars = "Age", .funs = gsub, pattern = "Y", replacement = "_yrs") %>%
  mutate_at(.vars = "Region/Area", .funs = gsub, pattern = "\\/", replacement = "-") %>% 
  mutate(colname = paste(AllenInstituteID, Age, `Region/Area`, sep = "_")) %>% 
  dplyr::select(-c(Agerange, Age, Description))

## Join the two files 
md = md %>% 
  left_join(.,md.excel, by = "colname", keep = TRUE)


## Process sheet 2
md.excel = read_excel(file.path(dir,"BrainSpan-additional.xlsx" ), sheet = 1, col_names = TRUE) %>% 
  as.data.frame() %>% 
  mutate_at(.vars="Internal ID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>% 
  dplyr::rename("Braincode" = "External ID")

md = merge(md, md.excel, by ="Braincode") 
md %<>% 
  dplyr::select(-c("Age.y", "colname.y", "colname.x", "Gender", AllenInstituteID,
                   "Region/Area", "age_for_mRIN", "sample name", "Internal ID")) %>% 
  dplyr::rename("Ethnicity"="Ethn.")  


dim(md)

table(md$SampleID)


message("Now adding Sample Information on expression matrix")

md = md[!duplicated(md[,c('column_num')]),]
md %>% 
  dplyr::arrange(column_num)

colnames(exp) <- md$SampleID
rownames(exp) <- rows.metadata$ensembl_gene_id
exp %<>% rownames_to_column("EnsemblID")




# Preprocessing of GTEx files 

## In the GTEx dir, you must ensure you have the following files 
## GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
## GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
## GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

gtex_dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/")

process_gtex(gtex_dir, outdir) # this function will save the files in the outdir 


### Extra annotations for bspan 

#### BrainSpan mRIN 

source("functions.R")


bspan.annot = read.csv("../../Results/Metadata/BrainSpan-metadata.csv", header=TRUE, check.names = FALSE, row.names =1) %>% 
  mutate(colname = paste(DonorName, Age, StructureAcronym, sep = "_"))


bspan.excel = read_excel("../../Annotations_and_updates/BrainSpan-additional.xlsx", sheet =2, col_names = TRUE, skip =1) %>% 
  as.data.frame() %>% mutate_at(.vars = "AllenInstituteID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>%
  mutate_at(.vars = "Age", .funs = gsub, pattern = "PCW", replacement = "_pcw") %>% 
  mutate_at(.vars = "Age", .funs = gsub, pattern = "M", replacement = "_mos") %>%
  mutate_at(.vars = "Age", .funs = gsub, pattern = "Y", replacement = "_yrs") %>%
  mutate_at(.vars = "Region/Area", .funs = gsub, pattern = "\\/", replacement = "-") %>% 
  mutate(colname = paste(AllenInstituteID, Age, `Region/Area`, sep = "_")) %>% 
  dplyr::select(-c(Agerange, Age, Description))



joined = bspan.annot %>% 
  left_join(.,bspan.excel, by = "colname", keep = TRUE)

bspan.excel.2 = read_excel("../../Annotations_and_updates/BrainSpan-additional.xlsx", sheet =1, col_names = TRUE) %>% 
  as.data.frame() %>% mutate_at(.vars="Internal ID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>% 
  rename("Braincode" = "External ID") 


joined = merge(joined, bspan.excel.2, by = "Braincode")
joined %<>% dplyr::select(-c("Age.y", "colname.y", "colname.x", "Gender", AllenInstituteID,
                             "Region/Area", "age_for_mRIN", "sample name", "Internal ID")) %>% 
  rename("Ethnicity"="Ethn." )



write.csv(joined, "../../Annotations_and_updates/BrainSpan-metadata-updated.csv")



# BrainSeq 


dir = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted")



clean_and_format(dir,"BrainSeq", outdir)




dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted")



clean_and_format(dir,"GTEx", outdir)

gtex = read.csv("../../../BrainData/Bulk/GTEx/Formatted/GTEx-metadata.csv")



dir = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted")
clean_and_format(dir, "PsychEncode", outdir)

