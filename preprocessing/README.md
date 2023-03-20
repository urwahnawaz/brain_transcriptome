# Pre-processing of RNA-seq datasets 

In order to preprocess the raw RNA-seq datafiles: 

On terminal:

```
git clone https://github.com/unawaz1996/brain_transcriptome.git

cd brain_transcriptome/preprocessing 
```

Open up Rstudio, or code editor of choice and ensure that you are working from the preprocessing directory 

In Rstudio, open file 

```
clean_and_format.R
```

This file sources to several functions that are written in the functions.R 


Just add the location of directories to where 
1) Raw data is stored 
2) Output directory 
3) Bulk dataset that's being formatted 

## BrainSpan

To preprocess BrainSpan data, you must ensure that the BrainSpan Dir contains:

* columns_metadata.csv
* rows_metadata.csv
* expression_matrix.csv

To retrieve these files: 

1) Go to BrainSpan Developmental Atlas 

2) Download the RNA-Seq Gencode v10 summarized to genes 

Two additional information files that contain information on the BrainSpan RNA-seq data were also retrieved for the pre-processing steps. These files are located in the `annotations` and will be automatically fetched once you run the script.  

```{r}
bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/")

clean_and_format(bspandir,"BrainSpan", outdir)
```

## GTEx 

To preprocess GTEx brain data - the following files must be in the GTEx input directory

* GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz 
* GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt 
* GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt



