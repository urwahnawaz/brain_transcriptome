source("libraries.R")



## Establish a set of cell-types 

ct.df = data.frame("Major_CT" = c(rep("Glia", 4), rep("PN", 4), rep("IN_CGE", 3), rep("IN_MGE",3)),
                "Subtype" = c("Astro", "Micro", "OPC", "Oligo", 
                              "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", 
                              "VIP", "ID2", "LAMP5_NOS1", 
                              "PV", "PV_SCUBE3", "SST"))


directory = file.path("/Users/urwah/Documents/PhD/Brain_transcriptome/Data/SN-data/major-dev-traj")
cell_types = list()
pattern = "/Users/urwah/Documents/PhD/Brain_transcriptome/Data/SN-data/major-dev-traj/major-dev-traj_"
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
signatures %<>% dplyr::select(-contains("Poor"))
    

## Signatures as they are

cpm_all =  apply(signatures, 2, function(x) {
    lib.size <- 10^6 / sum(x)
    x <- x * lib.size
    return(x)
})

cpm <- as.data.frame(cpm_all)
cpm

## rpkm 
## defined as 
### RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )

# get gene lengths 

## make TxDb from GTF file 
txdb <- makeTxDbFromGFF('/Users/urwah/Documents/PhD/Brain_transcriptome/Data/annotations/gencode.v19.annotation.gtf.gz')

## get gene information
all.genes <- genes(txdb)

## import your list of gene names

genes = signatures %>% rownames_to_column("genes") %>% 
    dplyr::select(genes) %>% 
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="")




ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"))
annotations <- dplyr::transmute(annotations, ensembl_gene_id, external_gene_name, gene_length = end_position - start_position)

# Filter and re-order gene.annotations to match the order in your input genes list
final.genes <- annotations %>% dplyr::filter(ensembl_gene_id %in% genes$genes)
final.genes <- final.genes[order(match(final.genes$ensembl_gene_id, genes$genes)),]; rownames(final.genes) <-NULL
dim(final.genes)
dim(signatures)



## get the length of each of those genes
my.genes.lengths <- width(all.genes[genes])
## put the names back on the lengths
names(my.genes.lengths) <- genes

## print lengths
print(my.genes.lengths)
signatures 

signatures %<>% 
    rownames_to_column("genes") %>%
    mutate_at(.vars = "genes", .funs = gsub, pattern = "\\--.*", replacement ="") %>% 
    column_to_rownames("genes")
    


signatures_rpkm = signatures[rownames(signatures) %in% final.genes$ensembl_gene_id,]
dim(signatures_rpkm)
expression.rpkm <- data.frame(sapply(signatures_rpkm, function(column) 10^9 * column / final.genes$gene_length / sum(column)))


rownames(expression.rpkm) = rownames(signatures_rpkm)
expression.rpkm

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


 
 
 
 for (files in directory){
     ct <- list.files(files, full.names = TRUE, pattern = "\\pseudo-bulk-cts_min10.csv$")
     
     for (j in ct){
         print(j)
         name = gsub(pattern, "", j)
         print(name)
     }}
 



