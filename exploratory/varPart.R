libs = c("dplyr", "ggplot2", "variancePartition", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


# load data 

bspan.exp = read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/BrainSpan/BrainSpan-exp.csv", row.names = 1) %>% column_to_rownames("EnsemblID")
colnames(bspan.exp) = gsub("X", "", colnames(bspan.exp))
colnames(bspan.exp) = gsub("\\.", "-", colnames(bspan.exp))
bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/BrainSpan/BrainSpan-metadata.csv")

form.bspan =  ~ AgeNumeric + (1|StructureAcronym) + (1|Sex) + (1|Period) + mRIN

perform.varPart = function(exp, md, form, region){
  
  # Filter expression
  exp = exp[apply(exp >= 1, 1, sum) >= 0.1*ncol(exp),]
  
  # Filter region
  md = md %>% dplyr::filter(Regions == region)
  exp = exp %>% dplyr::select(md$SampleID) 
  
  ## Run varpart 
  varpart = fitExtractVarPartModel(exp, form, md)
  
}

regions = unique(bspan.md$Regions)


for (r in regions){
  
  print(paste("Now performing variance partition for ", r))
  results = perform.varPart(bspan.exp, bspan.md, form.bspan, r)

  write.csv(file = paste0("BrainSpan-vatPart-", r,".csv"), results)
  
  
  
}


# BrainSeq

bseq.exp <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/BrainSeq/BrainSeq-exp.csv") %>% column_to_rownames("EnsemblID")
bseq.md <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/BrainSeq/BrainSeq-metadata.csv")

form.bseq <- ~ AgeNumeric + (1|Sex) + RIN +  (1|Diagnosis) + mito_Rate + rRNA_rate + totalAssignedGene + totalMapped
colnames(bseq.md)
regions = unique(bseq.md$Regions)


for (r in regions){
  
  print(paste("Now performing variance partition for ", r))
  results = perform.varPart(bseq.exp, bseq.md, form.bseq, r)
  
  write.csv(file = paste0("BrainSeq-vatPart-", r,".csv"), results)
  
}

## PsychEncode 


pe.exp <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/PsychEncode/PsychEncode-exp.csv", row.names =1 )  %>% column_to_rownames("EnsemblID")
rownames(pe.exp) = gsub("\\.[0-9]*$","",rownames(pe.exp))


#%>% column_to_rownames("EnsemblID")
pe.md <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/PsychEncode/PsychEncode-metadata.csv")

form.pe = ~ AgeNumeric + (1|Sex) +  (1|Diagnosis) + (1|Period) + (1|AgeInterval)

pe.exp = pe.exp[apply(pe.exp >= 1, 1, sum) >= 0.1*ncol(pe.exp),]

results = fitExtractVarPartModel(pe.exp, form.pe, pe.md)

 

## Gtex 

gtex.exp <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/Gtex/GTEx-exp.csv", row.names = 1)

gtex.md <- read.csv("/home/neuro/Documents/BrainData/Bulk/FormattedData/Gtex/GTEx-metadata.csv")
regions = unique(gtex.md$Regions)
regions

gtex.form <- ~ TotalNReads + rRNA_rate + (1|TypeofBatch) + (1|DateofBatch) + (1|BSS_Collection_side_code) + (1|AgeInterval) + (1|Sex) + (1|StructureAcronym) + IntergenicRate + RIN


for (r in regions){
  
  print(paste("Now performing variance partition for ", r))
  results = perform.varPart(gtex.exp, gtex.md, gtex.form, r)
  
  write.csv(file = paste0("GTEx-vatPart-", r,".csv"), results)
  
  
  
}