source("pre-processing.R")


# Preprocessing of BrainSpan files 

## In the bspan directory, you must ensure you have the following files 
## columns_metadata.csv
## rows_metadata.csv
## expression_matrix.csv

bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted")

process_bspan(bspandir, outdir)


# Preprocessing of GTEx files 

## In the GTEx dir, you must ensure you have the following files 
## GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
## GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
## GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

gtex_dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")

outDIR = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted")

process_gtex(gtex_dir, outDIR)