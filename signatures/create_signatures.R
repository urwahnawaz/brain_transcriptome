source("libraries.R")
source("functions.R")

## Establish a set of cell-types 

## type of signatures needed: 

## - All: 14 cell-types - Each dev -rpkm and cpm - DONE 
## - Adult and fetal - 14 cell-types - rpkm and cpm -DONE 

## - All -6 cell-types for each dev - inhibitory and excitatory summed - rpkm and cpm - DONE 
## - Adult and fetal -6 cell-types for each dev - inhibitory and excitatory summed - rpkm and cpm 



## All - 5 cell types - Each dev - all Neuronal summed - rpkm and cpm 
## Adult and fetal - 5 cell types - all neuronal summed- rpkm and cpm 

ct.df = data.frame("Major_CT" = c(rep("Glia", 4), rep("PN", 4), rep("IN_CGE", 3), rep("IN_MGE",3)),
                "Subtype" = c("Astro", "Micro", "OPC", "Oligo", 
                              "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", 
                              "VIP", "ID2", "LAMP5_NOS1", 
                              "PV", "PV_SCUBE3", "SST"), 
                "Major_CT_all" = c(rep("Glia", 4), rep("Neuron", 10)) )

dev.dt = data.frame("Period" = c("Prenatal", rep("Postnatal", 5)), 
                    "Stage" = c("Fetal", "Neonatal","Infancy", "Childhood", "Adolescence", "Adult"))


## Get data and create matrices 
directory = file.path("/home/neuro/Documents/BrainData/single-cell/herring/major-dev-traj")
cell_types = list()
pattern = "/home/neuro/Documents/BrainData/single-cell/herring/major-dev-traj/major-dev-traj_"
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


signatures = do.call(cbind, signatures)
signatures %<>%
    dplyr::select(-contains(c("Poor", "Vas"))) %>%
    rownames_to_column("genes") %>% 
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="") %>% 
    column_to_rownames("genes") %>% 
    dplyr::select(-Micro.Neonatal)

pan_neuronal_signatures = sapply(unique(ct.df$Major_CT)[-1], function(y){
    subtype <- ct.df %>%
        dplyr::filter(Major_CT == y) %>%
        dplyr::pull(Subtype)
    sapply(colnames(cell_types$Astro), function(x){
        cell_types %>%
            .[names(.) %in% subtype] %>%
            lapply(function(z){z[, colnames(z) == x]}) %>%
            do.call(cbind,.)  %>%
            rowSums() 
    }, simplify = FALSE)  %>%
        do.call(cbind,.) %>%
        set_colnames(paste(y, colnames(.), sep = ".")) 
}, simplify = FALSE) %>% do.call(cbind, .) %>% 
    as.data.frame() %>%
    rownames_to_column("genes") %>%
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="")

pan_neuronal_signatures = signatures %>% 
    dplyr::select(contains(c("Astro", "Micro", "Oligo", "OPC"))) %>% 
    rownames_to_column("genes") %>%
    full_join(pan_neuronal_signatures) %>% 
    column_to_rownames("genes")


all_neuro = sapply(unique(ct.df$Major_CT_all)[-1], function(y){
    subtype <- ct.df %>%
        dplyr::filter(Major_CT_all == y) %>%
        dplyr::pull(Subtype)
    sapply(colnames(cell_types$Astro), function(x){
        cell_types %>%
            .[names(.) %in% subtype] %>%
            lapply(function(z){z[, colnames(z) == x]}) %>%
            do.call(cbind,.)  %>%
            rowSums() 
    }, simplify = FALSE)  %>%
        do.call(cbind,.) %>%
        set_colnames(paste(y, colnames(.), sep = ".")) 
}, simplify = FALSE) %>% 
    do.call(cbind, .) %>% 
    as.data.frame() %>%
    rownames_to_column("genes") %>%
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="")


all_neuro = signatures %>% 
    dplyr::select(contains(c("Astro", "Micro", "Oligo", "OPC"))) %>% 
    rownames_to_column("genes") %>%
    full_join(all_neuro) %>% 
    column_to_rownames("genes")


## Normalization and filtering 
exp_thresh = 1 
pfc_signatures = list()


## CPM 

# 1) All CT and stages 
pfc_signatures$cpm_sig = cpm(signatures) %>% 
    .[which(apply(., 1, max) > exp_thresh),] 

# 2) All CT only Adult and fetal stages 
pfc_signatures$cpm_sig_adult_fetal = signatures %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    cpm(.) %>%
    .[which(apply(., 1, max) > exp_thresh),] 

# 3) pan neuronal - all stages 
pfc_signatures$cpm_panNeuro = cpm(pan_neuronal_signatures) %>% 
    .[which(apply(., 1, max) > exp_thresh),] 

# 4) pan neuronal -- adult and fetal stages 
pfc_signatures$cpm_panNeuro_adult_fetal = pan_neuronal_signatures %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    cpm(.) %>%
    .[which(apply(., 1, max) > exp_thresh),] 

# 5) all neuronal --all stages 
pfc_signatures$cpm_all_neuro = cpm(all_neuro) %>% 
    .[which(apply(., 1, max) > exp_thresh),]

# 6) all neuronal -- fetal and adult 
pfc_signatures$cpm_all_neuro_adult_fetal = all_neuro %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    cpm(.) %>%
    .[which(apply(., 1, max) > exp_thresh),] 


### RPKM 

length = getLength(rownames(signatures))

# 1)  All CT and stages 
pfc_signatures$rpkm_all = signatures %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]


# 2) All CT only Adult and fetal stages 
pfc_signatures$rpkm_fetal_adult = signatures %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]

# 3) pan neuronal - all stages
pfc_signatures$rpkm_panNeuro = pan_neuronal_signatures %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]

# 4) pan neuronal -- adult and fetal stages
pfc_signatures$rpkm_panNeuro_fetal_adult = pan_neuronal_signatures %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]

# 5) all neuro - all stages 
pfc_signatures$rpkm_all_neuro = all_neuro %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]

# 6) all neuro - fetal and adult 
pfc_signatures$rpkm_all_neuro_fetal_adult = all_neuro %>% 
    .[rownames(.) %in% length$ensembl_gene_id,] %>%
    dplyr::select(contains(c("Fetal", "Adult"))) %>%
    rpkm(., length$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]



### final signatures 

unique(colnames(signatures))
colnames(all_neuro)
ct.df
save(pfc_signatures, file= "../../Results/signatures/pfc_signatures.Rda")



pfc_signatures$rpkm_all_neuro
### maybe after the pyschencode 

