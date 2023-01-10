source("libraries.R")

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
                              "PV", "PV_SCUBE3", "SST"))


directory = file.path("/home/neuro/Documents/BrainData/single-cell/herring/major-dev-traj")
cell_types = list()
pattern = "/home/neuro/Documents/BrainData/single-cell/herring/major-dev-traj/major-dev-traj_"
signatures =  list()

exp_thresh <- 1
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
signatures %<>% dplyr::select(-contains(c("Poor", "Vas")))
    
signatures %<>% rownames_to_column("genes") %>% 
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="") %>% 
    column_to_rownames("genes") %>% 
    dplyr::select(-Micro.Neonatal)


pfc_signatures = list()


cpm =  apply(signatures, 2, function(x) {
    lib.size <- 10^6 / sum(x)
    x <- x * lib.size
    return(x)
})

cpm = cpm[which(apply(cpm, 1, max) > exp_thresh),]

pfc_signatures$cpm_all = cpm

## rpkm 

# Filter and re-order gene.annotations to match the order in your input genes list
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
annotations <- biomaRt::getBM(mart = ensembl, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"))
annotations <- dplyr::transmute(annotations, ensembl_gene_id, external_gene_name, gene_length = end_position - start_position)
annotations
final.genes <- annotations %>% dplyr::filter(annotations$ensembl_gene_id %in% rownames(signatures))
final.genes <- final.genes[order(match(final.genes$ensembl_gene_id, rownames(signatures))),]; rownames(final.genes) <-NULL

exp.rpkm = signatures[rownames(signatures) %in% final.genes$ensembl_gene_id,]
expression.rpkm <- data.frame(sapply(exp.rpkm, function(column) 10^9 * column / final.genes$gene_length / sum(column)))
rownames(expression.rpkm) = rownames(exp.rpkm)
expression.rpkm = expression.rpkm[which(apply(expression.rpkm, 1, max) > exp_thresh),]
pfc_signatures$rpkm_all = expression.rpkm



## Fetal and adult 

fetal_cpm = signatures %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>% apply(., 2, function(x) {
    lib.size <- 10^6 / sum(x)
    x <- x * lib.size
    return(x)
})

fetal_cpm = fetal_cpm[which(apply(fetal_cpm, 1, max) > exp_thresh),]
pfc_signatures$all_adult_fetal_cpm = fetal_cpm


rpkm_fetal = exp.rpkm %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>% 
    data.frame(sapply(., function(column) 10^9 * column / final.genes$gene_length / sum(column)))
rpkm_fetal = rpkm_fetal[which(apply(rpkm_fetal, 1, max) > exp_thresh),]

pfc_signatures$all_adult_fetal_rpkm = rpkm_fetal

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


summed = sapply(unique(ct.df$Major_CT)[-1], function(y){
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
        set_colnames(paste(y, colnames(.), sep = "_")) 
}, simplify = FALSE)

# as long as number of column number and row number matches 
pan_neuronal_signatures = do.call("cbind", summed) %>% 
    as.data.frame() %>%
    rownames_to_column("genes") %>%
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="")


head(pan_neuronal_signatures,10)

pan_neuronal_signatures = signatures %>% 
    dplyr::select(contains(c("Astro", "Micro", "Oligo", "OPC"))) %>% 
    rownames_to_column("genes") %>%
    full_join(pan_neuronal_signatures) %>% column_to_rownames("genes")


pan_neuronal_signatures
PN_cpm = pan_neuronal_signatures %>% 
    apply(., 2, function(x) {
        lib.size <- 10^6 / sum(x)
        x <- x * lib.size
        return(x)
    })

PN_cpm = PN_cpm[which(apply(PN_cpm, 1, max) > exp_thresh),]

pfc_signatures$PN_all_cpm = PN_cpm

## RPKM 
exp.rpkm.PN = pan_neuronal_signatures[rownames(pan_neuronal_signatures) %in% final.genes$ensembl_gene_id,]
exp.rpkm.PN

data.frame(sapply(exp.rpkm.PN, function(column) 10^9 * column / final.genes$gene_length / sum(column)))


exp.rpkm.PN = exp.rpkm.PN[which(apply(exp.rpkm.PN, 1, max) > exp_thresh),]
pfc_signatures$PN_all_rpkm = exp.rpkm.PN

exp.rpkm.PN

## PN - fetal and adult 

PN_fetal_adult = 
    pan_neuronal_signatures[rownames(pan_neuronal_signatures) %in% final.genes$ensembl_gene_id,]  %>% 
    dplyr::select(contains(c("Fetal", "Adult"))) %>% 
    data.frame(sapply(., function(column) 10^9 * column / final.genes$gene_length / sum(column)))
PN_fetal_adult

rpkm_fetal = rpkm_fetal[which(apply(rpkm_fetal, 1, max) > exp_thresh),]
rownames(rpkm_fetal) = rownames(rpkm_fetal)



y <- strsplit(rownames(trial), "--")
ensid <- sapply(y, "[", 1)
ensid




