## Deconvolve HDBR 

source("libraries.R")
source("functions.R")
load("../../DeconRNAShiny/sigsBrain.rda")


hdbr.md = read.csv("/home/neuro/Documents/BrainData/Bulk/HDBR/Formatted/HDBR-metadata.csv", header=TRUE, check.names = FALSE)
hdbr.exp  = read.csv("/home/neuro/Documents/BrainData/Bulk/HDBR/Formatted/HDBR-exp.csv", 
                     header=TRUE, check.names = FALSE, row.names = 1) %>% 
    column_to_rownames("EnsemblID")


## Run dtangle using MultiBrain 
hbdr.decon = run_dtg(hdbr.exp ,sigsBrain$MB) #%>%
# as.data.frame()


gof_res = write.gof(hdbr.exp, hbdr.decon, 
                    signatureUsed = sigsBrain$MB)


gof_res %>% 
    mutate(col = "HDBR") %>%
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% 
 
    left_join(hdbr.md, by = "SampleID") %>%
    drop_na(Regions) %>%
    ggplot(aes(Regions, r)) + geom_violin() 

head(gtex.decon)

### Brainseq 
bseq.md = read.csv("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-metadata.csv", header=TRUE, check.names = FALSE)
bseq.exp  = read.csv("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-exp.csv", 
                     header=TRUE, check.names = FALSE, row.names = 1) 


## Run dtangle using MultiBrain 
bseq.decon = run_dtg(bseq.exp ,sigsBrain$MB) #%>%
# as.data.frame()


gof_res = write.gof(bseq.exp, bseq.decon, 
                    signatureUsed = sigsBrain$MB)


gof_res %>% 
    mutate(col = "Bseq") %>%
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% 
    
    left_join(bseq.md, by = "SampleID") %>%
    drop_na(Regions) %>%
    ggplot(aes(Regions, r)) + geom_violin() +
    facet_grid(~Period)

bseq.decon %>% 
    rownames_to_column("SampleID") %>% 
    left_join(bseq.md, by = "SampleID") %>%
    dplyr::filter(Period == "Prenatal") %>% 
    ggplot(aes())
