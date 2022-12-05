# Pre-processing of RNA-seq datasets 

In order to preprocess the raw RNA-seq datafiles: 

On terminal:

```
git clone https://github.com/unawaz1996/brain_transcriptome.git

cd brain_transcriptome/preprocessing 
```

Open up Rstudio, or code editor of choice and ensure that you are working from the preprocessing direcotry 

In Rstudio, open file 

```
clean_and_format.R
```

This file sources to several functions that are written in the pre-processing.R script. 

Just add the location of directories of where the raw data is stored 

# BrainSpan 

To preprocess the BrainSpan data, you must ensure your directory of raw data contains 

* columns_metadata.csv
* rows_metadata.csv
* expression_matrix.csv

```
source("pre-processing.R")

bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted")

process_bspan(bspandir, outdir)
```

This script will then save the BrainSpan-exp.cvs and BrainSpan-metadata.csv file into the results directory you provided. 


# GTEx

To preprocess the GTEx files, you must ensure you have the following files in your directory: 
* GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

* GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

* GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

```
source("pre-processing.R")

gtex_dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")

outDIR = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted")

process_gtex(gtex_dir, outDIR)

```


