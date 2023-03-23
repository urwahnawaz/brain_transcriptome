# Pre-processing of RNA-seq datasets

In order to preprocess the raw RNA-seq datafiles: 

On terminal:

```
git clone https://github.com/unawaz1996/brain_transcriptome.git

cd brain_transcriptome/preprocessing 
```

## Bulk RNA-seq datasets 

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


### BrainSeq 
To preprocess BrainSeq data, the following two files are needed:

BrainSeq2 [expression and metadata file](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata)
These files are R objects. To run the preprocessing steps, change file path to where the R objects are located

```
# Change to your own input directory where BrainSeq files are located 
dir = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq")

# Change to your own out directory
outdir = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted")

clean_and_format(dir,"BrainSeq", outdir)

```

### BrainSpan

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

### GTEx 

To preprocess GTEx brain data - the following files must be in the GTEx input directory

* GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz 
* GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt 
* GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

```
# Change to your own input directory where GTEx files are located 
dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")

# Change to your own out directory
outdir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted")

clean_and_format(dir,"GTEx", outdir)

```

### PsychEncode 

Retriving and preprocessing files from PsychEncode will require access to the Synapse Portal for the metadata files. The expression matrix can be downloaded from [here](http://resource.psychencode.org/Datasets/Derived/DER-02_PEC_Gene_expression_matrix_TPM.txt)

To preprocess the data - the following files must be in the input directory: 

* [DER-02_PEC_Gene_expression_matrix_TPM.txt](http://resource.psychencode.org)

* [DER-24_Cell_fractions_Normalized.xlsx](http://resource.psychencode.org/Datasets/Derived/SC_Decomp/DER-24_Cell_fractions_Normalized.xlsx)

* Metadata file retrieved from *Synapse* - Ensure that the metadata file starts with `Job`


```
# Change to your own input directory where PsychEncode files are located 

dir = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode")

# Change to your own out directory
outdir = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted")

clean_and_format(dir,"PsychEncode", outdir)

```




## Single-cell data pre-processing 

In order to preprocess the single-nuclues RNA-seq datasets, open the `Script.R` in Rstudio and change all input and output file locations to relevant local directories and run the analysis. 


### Velmeshev et al 

To process data from Velmeshev et al (2019), the data must be acquired through the [UCSC cell browser](https://cells.ucsc.edu/?ds=autism). The following files must be in the input directory:

* barcodes.tsv
* genes.tsv
* matrix.mtx 
* meta.txt

These can be downloaded through [here](https://cells.ucsc.edu/autism/rawMatrix.zip) 

```
## Change the following to your input directory
vel = file.path("/home/neuro/Documents/BrainData/single-cell/velmeshev/Velmeshev2019/")

## Change the following to your output directory
vel.outfile = file.path("/home/neuro/Documents/BrainData/single-cell/velmeshev/FormattedData")
```

# Subsetting of data and determining drivers of variation 

The `md_correlation.html` file goes through correlation plots of the metadata variables from bulk RNA-seq datasets and goes through performing the `variancePartition` analysis. Please refer to this document for more information.
