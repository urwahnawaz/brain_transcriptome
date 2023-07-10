## Deconvolve GTEx data 

source("libraries.R")
source("functions.R")
load("../../DeconRNAShiny/sigsBrain.rda")


gtex.md = read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-metadata-subset.csv", header=TRUE, check.names = FALSE)
gtex.exp  = read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-exp.csv", header=TRUE, check.names = FALSE, row.names = 1)


## Run dtangle using MultiBrain 
gtex.decon = run_dtg(gtex.exp ,sigsBrain$MB) #%>%
   # as.data.frame()


gof_res = write.gof(gtex.exp, gtex.decon, 
                    signatureUsed = sigsBrain$MB)


gof_res %>% 
    mutate(col = "GTEx") %>%
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% 
    mutate(SampleID = gsub("\\.", "-", SampleID)) %>% 
    left_join(gtex.md, by = "SampleID") %>%
    ggplot(aes(Regions, r)) + geom_violin() 

head(gtex.decon)


write.csv(gtex.decon, file = "/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-decon-MB.csv")
