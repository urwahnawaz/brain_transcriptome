library(edgeR)
library(dplyr)
library(data.table)
library(tidyverse)
library(magrittr)
library(plyr)
library(biomaRt)
library(EDASeq)
library(stringr)

# sum columns 


## Establish a set of cell-types 

ct.df = data.frame("Major_CT" = c(rep("Glia", 4), rep("PN", 4), rep("IN_CGE", 3), rep("IN_MGE",3)),
                "Subtype" = c("Astro", "Micro", "OPC", "Oligo", 
                              "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", 
                              "VIP", "ID2", "LAMP5_NOS1", 
                              "PV", "PV_SCUBE3", "SST"))


directory = file.path("/Volumes/share/mnt/Data0/PROJECTS/GWAS_Enrichment/GAVIN/Data/Expression/major-dev-traj")
cell_types = list()
pattern = "/Volumes/share/mnt/Data0/PROJECTS/GWAS_Enrichment/GAVIN/Data/Expression/major-dev-traj/major-dev-traj_"
signatures =  list()


for (files in directory){
    ct <- list.files(files, full.names = TRUE, pattern = "\\pseudo-bulk-cts_min10.csv$")
    
    for (j in ct){
        
        ct_file = read.csv(j, header= TRUE, row.names = 1)
        prefixes = unique(sub("\\..*", "", colnames(ct_file)))
        name = gsub(pattern, "", j)
        name = gsub("_pseudo-bulk-cts_min10.csv", "", name)
        message("Now manipulating data for ", paste0(name))
        
        for (prefix in prefixes){
            column = ct_file %>% dplyr::select(contains(prefix)) 
            if (ncol(column) > 1) {
                message("All good, keeping columns for ", paste0(prefix))
                single_col = NULL
            } else {
                message(paste0(prefix), " only contains a single column")
                #name = colnames(column)
                single_col = ct_file %>% as.data.frame() %>% dplyr::select(colnames(column))
                ct_file = ct_file %>% as.data.frame() %>%  dplyr::select(-colnames(column))
            }
            
        }
        
        message("Now summing for ", paste0(name))
        signature = sapply(prefixes, function(x)rowSums(ct_file[,startsWith(colnames(ct_file), x)]))
        
        
        
        if (is.null(single_col)) {
        } else {
            signature = data.frame(signature, single_col)
        }
        final_sig = signature 
        final_sig %>% set_colnames(paste(name, colnames(final_sig), sep = "_"))
         
        ## Add all different cell types 
        cell_types[[paste0(name)]] = signature 
        signatures[[paste0(name)]] = final_sig
     }
    
}



do.call(cbind, signatures)

## Signatures as they are
signatures =  list()
#signatures$all= 
    
    cell_types %>% 
        do.call(cbind,.) %>% 
        set_colnames(paste(names(cell_types), colnames(.), sep = "_")) 
        
    #do.call(cbind,.)%>% 

    
lapply(cell_types, function(x){
    set_colnames(paste(names(x), colnames(x), sep = "_"))
})
colnames()
head(cell_types)
matrix = do.call("cbind", cell_types)
matrix


ct.df
summed = sapply(unique(ct.df$Major_CT)[-1], function(y){
    subtype <- ct.df %>%
        dplyr::filter(Major_CT == y) %>%
        dplyr::pull(Subtype)
    sapply(colnames(cell_types$Astro), function(x){
        cell_types %>%
            .[names(.) %in% subtype] %>%
            lapply(function(z){z[, colnames(z) == x]}) %>%
            do.call(cbind,.)  %>%
            rowSums() %>%
            head()
    }, simplify = FALSE)  %>%
        do.call(cbind,.) %>%
        set_colnames(paste(y, colnames(.), sep = "_")) 
}, simplify = FALSE)

# as long as number of column number and row number matches 
summed

## Example 



## creating cpm signatures 
#### sigs is cpm based 
sigs = list()
all = cpm(cell_types_signature)
dim(all)
filt = rowSums(all >= 1 ) >= 6
all = all[filt,]
sigs$all_CT = all


## 14 x all Fetal and adult only 
sigs_fetal_adult = cell_types_signature %>% as.data.frame() %>% 
    dplyr::select(contains("Fetal"), contains("Adult"))

sigs_fetal_adult = cpm(sigs_fetal_adult)
filt = rowSums(sigs_fetal_adult >= 1 ) >= 2
sigs_fetal_adult= sigs_fetal_adult[filt,]
dim(sigs_fetal_adult)
sigs$fetal_adult = sigs_fetal_adult




## Sum all neuronal cell-types per category 





### Downloadable from 10x genomics 
# Need to get gene length

# convert to rpkm 
# first round filer = > 1 rpkm per cell type 


#second round filer = > 1 rpkm per stage 


y <- strsplit(rownames(trial), "--")
ensid <- sapply(y, "[", 1)
ensid


human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
 gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=ensembl_list, mart=human)
 gene_coords$size=gene_coords$end_position - gene_coords$start_position
 gene_coords
 
 
 
 



