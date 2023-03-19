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

For example, to preprocess BrainSpan data, you must ensure that the BrainSpan Dir contains:

* columns_metadata.csv
* rows_metadata.csv
* expression_matrix.csv

```{r}
bspandir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")

outdir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/")

clean_and_format(bspandir,"BrainSpan", outdir)
```



