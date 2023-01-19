suppressPackageStartupMessages({
  library(magrittr)
  library(Matrix)
  library(data.table)
  library(gdata)
  library(DT)
  library(tidyverse)
  library(edgeR)
})

## -- load data 
path = file.path("/home/neuro/Documents/BrainData/single-cell/velmeshev")
meta <- read.table(file.path(path, "meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)

counts <- readMM(file.path(path, "matrix.mtx")) 
genes <- read_tsv(file.path(path, "genes.tsv"), col_names = FALSE)
barcodes <- read_tsv(file.path(path, "barcodes.tsv"), col_names=FALSE)


colnames(counts) = barcodes$X1 ### cells as colnames 
rownames(counts) = genes$X1 ### genes for rownames 


counts %<>%  ## converting dgTMatrix option into a dataframe 
  as.matrix() %>%
  as.data.frame() %>% 
  cpm()
dim(counts)[1:10]



autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")



ASD_list <- unique(autism$individual)


zscore <- data.frame(matrix())

for (i in ASD_list) {
  ### subsetting for the first individual
  ### get all the information out from the meta data 
  ## by using dplyr we are retriving all cells that belong to an ASD individual
  ## and controls 
  name <- i
  
  message("subsetting from meta data for"," ", i)
  
  meta_ASD <- meta %>% 
    rownames_to_column("sample_names") %>%
    dplyr::filter(individual == name | diagnosis == "Control") %>% 
    column_to_rownames("sample_names")    
  
  ### subset the counts based on the ASD meta data 
  
  message("Subsetting counts for"," ",i)
  counts_ASD <- counts %>% 
    as.data.frame() %>% 
    dplyr::select(as.character(rownames(meta_ASD)))
  
  idList <- meta_ASD %>% 
    rownames_to_column("cellID") %>%
    mutate(f = paste(cluster, individual, sep = "_")) %>% 
    split(f = .$f) %>% lapply(extract2, "cellID")    ### by using lapply we sum the counts 
  
  message("Now calculating averages based on individual and celltype for"," ", i)
  
  counts_averaged <- lapply(idList, function(x){rowMeans(counts[,x])}) 
  counts_averaged <- as.data.frame(do.call(rbind, counts_averaged))
  counts_averaged <- as.data.frame(t(counts_averaged[,-1]))    ## retrieve cell-type information from this
  
  patterns <- unique(sub("_.*", "", colnames(counts_averaged)))    
  
  message("Calculating zscore for each celltype for"," ", i)
  for (m in patterns){
    cell_type <- counts_averaged %>% 
      dplyr::select(matches(m))
    zscore_calc_cellType <- as.data.frame(t(scale(cell_type)))
    
    zscore_results <- as.data.frame(t(zscore_calc_cellType)) %>%
      dplyr::select(matches(as.character(i)))
    message("zscore calculate for cluster", " ", m)
    zscore <- cbind(zscore, zscore_results)
  }
  
  message("All finished for"," ", i, ":)")
}




zscore_ASD = function(ASD_list, counts, meta) {
  
  ### Loading all the necessary libraries for this analysis 
  suppressPackageStartupMessages({
    library(magrittr)
    library(Matrix)
    library(data.table)
    library(gdata)
    library(DT)
    library(tidyverse)
    library(edgeR)
  })
  
  counts <- counts 
  counts %<>%  ## converting dgTMatrix option into a dataframe 
    as.matrix() %>%
    as.data.frame() %>% 
    cpm()
  zscore <- data.frame(matrix())
  for (i in ASD_list) {
    ### subsetting for the first individual
    ### get all the information out from the meta data 
    ## by using dplyr we are retriving all cells that belong to an ASD individual
    ## and controls 
    name <- i
    message("subsetting from meta data for"," ", i)
    meta_ASD <- meta %>% 
      rownames_to_column("sample_names") %>%
      dplyr::filter(individual == name | diagnosis == "Control") %>% 
      column_to_rownames("sample_names")    
    
    ### subset the counts based on the ASD meta data 
    
    message("Subsetting CPM normalised counts for"," ",i)
    
    counts_ASD <- counts %>% 
      as.data.frame() %>% 
      dplyr::select(as.character(rownames(meta_ASD)))
    
    ### creating an Id list
    ## this id list contains informaion of individual and clusters
    ## both of which will be used to sum counts (i.e counts for cells that belong to a cluster 
    ## from a specific individual will be summed)
    
    idList <- meta_ASD %>% 
      rownames_to_column("cellID") %>%
      mutate(f = paste(cluster, individual, sep = "_")) %>% 
      split(f = .$f) %>% lapply(extract2, "cellID")    ### by using lapply we sum the counts 
    
    message("Now summing CPM counts based on individual and celltype for"," ", i)
    
    counts_summed <- lapply(idList, function(x){rowSums(counts[,x])}) 
    ## lapply returns a list, so to flatten the list we use the following function
    counts_summed <- as.data.frame(do.call(rbind, counts_summed)) 
    counts_summed <- as.data.frame(t(counts_summed[,-1]))    
    
    ## retrieve cell-type information from this
    patterns <- unique(sub("_.*", "", colnames(counts_summed)))    
    
    message("Calculating zscore for each celltype for"," ", i)
    ## we want to calculate zscores for each cell type for each ASD individual 
    ## so we subset cells based if they belong to a specific cell-type
    
    for (m in patterns){
      cell_type <- counts_summed %>% 
        dplyr::select(matches(m))
      
      zscore_calc_cellType <- as.data.frame(t(scale(cell_type)))
      
      # we are just interested in the zscores for the ASD individual so we susbet that result 
      # and save into a DF
      zscore_results <- as.data.frame(t(zscore_calc_cellType)) %>%
        dplyr::select(matches(as.character(i)))
      message("zscore calculate for cluster", " ", m)
      zscore <- cbind(zscore, zscore_results)
    }
    
    message("All finished for"," ", i, ":)")
  }
}


head(zscore)
write.csv(zscore, "zscores_summed.csv")


## 
random_list = read.csv("/home/neuro/Documents/Archive/CP_transcriptome_analysis/Data/Gene_lists/CP_genes.csv", header= TRUE)
CP = zscore[rownames(zscore) %in% random_list$Gene, ]
dim(CP)
CP


CP  <- CP %>%
  rownames_to_column("geneID") %>%
  melt()


clusters <- unique(sub("_.*", "", CP$variable))
head(patterns)

CP %<>% 
  mutate(cluster = sub("_.*", "", CP$variable)) %>% 
  mutate(individual = sub(".*_", "", CP$variable))

table(CP$individual)

ggplot(CP, aes(x=cluster, y=log2(value), fill = cluster)) + geom_boxplot() + 
  ylab("zscore") + theme_bw() + theme(legend.position = "none", 
                                      axis.text.x = element_text(angle = 90))

