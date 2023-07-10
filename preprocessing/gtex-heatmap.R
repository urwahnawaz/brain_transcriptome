# Replicating GTEx heatmap for ZMYND8

library(CePa)
library(EnsDb.Hsapiens.v86)
library(recount3)
library(DT)
library(dplyr)
library(tidyr)
source("functions.R")
library(pander)
library(gridExtra)
library(variancePartition)
library(corrplot)
library(edgeR)
library(pheatmap)
library(viridis)

#gtex.transcript =  read.table("/home/neuro/Documents/BrainData/Bulk/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",  sep="\t", skip=2, header=TRUE)

head(gtex.transcript)[1:10]

#file= read.gct("/home/neuro/Documents/BrainData/Bulk/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")

gtex_trans = read.gct("/home/neuro/Documents/BrainData/Bulk/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")
colnames(gtex_trans ) <- gsub("\\.", "-", colnames(gtex_trans))
gtex.bH.md.subset = read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-metadata-subset.csv", header=TRUE, check.names = FALSE)
txdf = transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene = as.data.frame(txdf[,c("gene_id", "tx_id")])

gtex_trans_brain = gtex_trans %>% 
  as.data.frame() %>%
  dplyr::select(contains(gtex.bH.md.subset$SampleID))

zmynd8_heatmap = gtex_trans_brain %>% 
  rownames_to_column("tx_id") %>%
  mutate_at(.vars = "tx_id", .funs = gsub, pattern = "\\.[0-9]*$", replacement = "") %>% 
  left_join(tx2gene, by ="tx_id") %>% 
  dplyr::filter(gene_id == "ENSG00000101040") %>% 
  column_to_rownames("tx_id") %>% 
  dplyr::select(-gene_id)

gtex.bH.md.subset %<>% 
  dplyr::arrange(StructureAcronym)

sub_samp_ordered <- zmynd8_heatmap[,gtex.bH.md.subset$SampleID]


pheatmap(sub_samp_ordered, scale = "none", cluster_rows = FALSE, 
         cluster_cols= FALSE)


## Average by Region 

# sapply(unique(gtex.bH.md.subset$StructureAcronym), function(y){
#   structure <- gtex.bH.md.subset %>%
#     dplyr::filter(StructureAcronym == y) %>%
#     dplyr::pull(SampleID)})
#   sapply(colnames(sub_samp_ordered), function(x){
#     sub_samp_ordered %>%
#       .[names(.) %in% structure] %>%
#       lapply(function(z){z[, colnames(z) == x]}) %>%
#       do.call(cbind,.)  %>%
#       rowSums() 
#   }, simplify = FALSE)})


# #%>%
#     do.call(cbind,.) %>%
#     set_colnames(paste(y, colnames(.), sep = ".")) 
# }, simplify = FALSE) %>% do.call(cbind, .) %>% 
#   as.data.frame() 
# 
# #%>%
#   rownames_to_column("genes") %>%
#   mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="")
#   
  
#### Solution online - works!

                                                                                                                            
average = with(gtex.bH.md.subset, setNames(StructureAcronym, SampleID))[names(sub_samp_ordered)[col(sub_samp_ordered)]]

tapply(unlist(sub_samp_ordered), list(row(sub_samp_ordered), average ), mean)

a <- with(df_category, setNames(Category, Col_name))[names(df)[col(df)]]
  tapply(unlist(df), list(row(df), a), mean)
  
  
### Recreate ZMYND8 transcript heatmap from GTEx using all data 
  

dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/")
annot= read.csv(paste0("../annotations/","GTEx","-metadata-annot.csv"))
attributes = list.files(dir, full.names = TRUE, pattern = "\\SampleAttributesDS.txt") # Sample attributes contains sample level information
phenotype = list.files(dir, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt")
md = read_tsv(attributes, col_types = c('.default' = 'c')) %>% 
    #filter(SMTS == 'Brain') %>% 
    mutate(SUBJID = sapply(str_split(SAMPID, pattern = "-"), function(x) paste(x[1:2], collapse = '-'))) %>%
    left_join(read_tsv(phenotype, col_types = c('.default' = 'c')))  %>% as.data.frame()
colnames(md) = annot$BITColumnName[match(colnames(md), annot$OriginalMetadataColumnName)]
  
zmynd8_heatmap = gtex_trans %>%
  as.data.frame() %>%
  rownames_to_column("tx_id") %>%
  mutate_at(.vars = "tx_id", .funs = gsub, pattern = "\\.[0-9]*$", replacement = "") %>% 
  left_join(tx2gene, by ="tx_id") %>% 
  dplyr::filter(gene_id == "ENSG00000101040") %>% 
  column_to_rownames("tx_id") %>% 
  dplyr::select(-gene_id) %>% 
  dplyr::select(contains(md$SampleID))

dim(zmynd8_heatmap)

StructureAcronym


average_with_age = with(md, setNames(StructureAcronym, AgeInterval, SampleID))[names(zmynd8_heatmap)[col(zmynd8_heatmap)]]

average_exp = tapply(unlist(zmynd8_heatmap), list(row(zmynd8_heatmap), average ), mean)

rownames(average_exp) = rownames(zmynd8_heatmap)


pheatmap(log(average_exp + 0.05), scale = "none", cluster_rows = FALSE, 
         cluster_cols= FALSE,cellwidth = 20, cellheight = 20, border_color = "black", 
         color= magma(10))


a <- with(df_category, setNames(Category, Col_name))[names(df)[col(df)]]
tapply(unlist(df), list(row(df), a), mean)



  