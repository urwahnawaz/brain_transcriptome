## dtangle 
source("libraries.R")


annot = colnames(signatures) %>% 
    as.data.frame() %>% 
    set_colnames("Sample") %>% 
    mutate(Cell_type = gsub("\\..*", "", Sample))



pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(anno$CellType == all_cell_type[i])
})