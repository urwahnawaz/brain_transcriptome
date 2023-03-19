libs = c("dplyr", "ggplot2", "variancePartition", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


# load data 

bspan.exp = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/BrainSpan-exp.csv", row.names = 1, check.names = FALSE) %>%
    column_to_rownames("EnsemblID")


colnames(bspan.exp) = gsub("X", "", colnames(bspan.exp))
colnames(bspan.exp) = gsub("\\.", "-", colnames(bspan.exp))
bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/BrainSpan-metadata.csv")

form.bspan =  ~ AgeNumeric + (1|StructureAcronym) + (1|Sex) + (1|Period) + mRIN + RIN


## Run variance parition

perform.varPart = function(exp, md, form, region){
    
    exp = exp[apply(exp >= 1, 1, sum) >= 0.1*ncol(exp),]
    
    
    if(missing(region)){
        varpart = fitExtractVarPartModel(exp, form, md)
    } else {
        md = md %>% dplyr::filter(Regions == region)
        exp = exp %>% dplyr::select(md$SampleID)
        varpart = fitExtractVarPartModel(exp, form, md)
        
    }
 
  return(varpart)
}





bspan.res = perform.varPart(bspan.exp, bspan.md, form.bspan)



bspan.res + theme_bw()

bspan.res= sortCols(bspan.res)

plotVarPart(bspan.res)

varPart.plot= bspan.res %>% melt() %>% 
    ggplot(aes(variable, value)) +geom_violin(scale = "width", color = "white", width = 0.9, fill = alpha('#F2AF1D', 0.75)) +
    geom_boxplot(width=0.07,fill= "grey" , outlier.colour = "grey") +
    xlab("") + ylab("Variance explained (%)") + theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                                      plot.background = element_rect(fill = "transparent",colour = NA),
                                                      axis.text.x = element_text(colour = "white", face = "bold"), 
                                                      axis.text.y = element_text(color = "white", face = "bold"), 
                                                      axis.title.y = element_text(color = "white", face = "bold"))

varPart.plot + scale_y_reverse()

ggsave(file="../../Results/exploratory/varPart_exp.svg", plot=varPart.plot, width=8.3, height=6.25, units = "in")


plot.eg = bspan.res["ENSG00000007372",] %>% melt() %>% 
    ggplot(aes(x= reorder(variable, value), y = value)) + geom_bar(stat = "identity", color = "#F16160", fill = alpha("#F59392", 0.75), size =1) + 
    theme_bw() + ylab("Proportion of variance explained") + xlab("") + coord_flip() +
    theme(axis.text.y = element_text(hjust=1)) +
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour = NA))

plot.eg 
ggsave(file="../../Results/exploratory/varPart_exp_PAX6.svg", plot=plot.eg , width=8.3, height=6.25, units = "in")



bseq_varPart["ENSG00000116273",]


?ggsave
plot +
    theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA), 
    legend.position = "none") + theme_classic()

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