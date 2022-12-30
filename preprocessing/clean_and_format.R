source("pre-processing.R")


# Preprocessing of BrainSpan files 

## In the bspan directory, you must ensure you have the following files 
## columns_metadata.csv
## rows_metadata.csv
## expression_matrix.csv

bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/")

process_bspan(bspandir, outdir) # This function will save the files in the outdir


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