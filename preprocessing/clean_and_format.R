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

dim(joined)
colnames(joined)
pdf("../../Results/mRIN/BrainSpan_RIN_vs_mRIN.pdf", height = 6, width = 6)
joined %>% 
  ggscatter(x="RIN", y ="mRIN", cor.coef = TRUE, add = "reg.line", 
            conf.int = TRUE, add.params = list(color = "blue",
                                               fill = "lightgray")) + 
  theme_bw() + geom_hline(yintercept = -0.033, linetype="dashed") +
  geom_vline(xintercept = 7, linetype="dashed", color = "red") + ggtitle("BrainSpan RIN vs mRIN")
dev.off()

bspan.excel.2 = read_excel("../../Annotations_and_updates/BrainSpan-additional.xlsx", sheet =1, col_names = TRUE) %>% 
  as.data.frame() %>% mutate_at(.vars="Internal ID", .funs = gsub, pattern = "\\.", replacement = "\\_") %>% 
  rename("Braincode" = "External ID") 
dim(bspan.excel.2)




joined = merge(joined, bspan.excel.2, by = "Braincode")
joined %<>% dplyr::select(-c("Age.y", "colname.y", "colname.x", "Gender", AllenInstituteID,
                             "Region/Area", "age_for_mRIN", "sample name", "Internal ID")) %>% 
  rename("Ethnicity"="Ethn." )


dim(joined)

write.csv(joined, "../../Annotations_and_updates/BrainSpan-metadata-updated.csv")


### single nucleus data cleaning and preprocessing 
read.csv("../../Data/SN Formatted /Velmeshev-metadata.csv", header= TRUE) %>% 
  mutate(Age = paste0(AgeNumeric, " yrs")) %>% 
  mutate(AgeIntervals = add_feature(.$Age, age_intervals),
         Regions = add_feature(.$StructureAcronym, regions)) %>% 
  apply(.,2,as.character) %>% 
  write.csv(., "../../Data/SN Formatted /Velmeshev-metadata_updated.csv")

Vel = read.csv("../../Data/SN Formatted /Velmeshev-metadata_updated.csv", header=TRUE)
head(Vel)
table(Vel$AgeIntervals)

HCA =read.csv("../../Data/SN Formatted /HCA-metadata_fixed.csv", header= TRUE) 
table(HCA$region)
dim(HCA)
