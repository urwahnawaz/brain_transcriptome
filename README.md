# BITHub

This repository contains the code used for preprocessing transcriptomic files used for Brain Intergrative Transcriptome database (BITHub)


## Datasets 

### Bulk RNA-seq datasets 


| Dataset | Raw Files | Formatted file |
|---------|-----------|----------------|
| BrainSpan | columns_metadata.csv <br> rows_metadata.csv <br> expression_matrix.csv | BrainSpan-metadata.csv <br> BrainSpan-exp.csv |
| BrainSeq | rse_gene_unfiltered.Rdata <br> methprop_pd.Rdata | BrainSeq-metadata.csv <br> BrainSeq-exp |
| GTEx | GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt <br> GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt <br> GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | GTEx-metadata.csv <br> GTEx-exp.csv |
| PsychEncode | Job-154061626746631166818172835.csv <br> DER-24_Cell_fractions_Normalized.xlsx <br> DER-01_PEC_Gene_expression_matrix_TPM.txt | PsychEncode-metadata.csv <br> PsychEncode-exp.csv | 

### Single-cell RNA-seq datasets
|DataSet | Raw Files | Formatted files |
|--------|-----------|---------------- |
| Velmesshav et al |||
| Human Cell Atlas |||
#### Scripts 



### Exploratory analysis


### How to run to update website 

#### Changing colours
* git clone https://github.com/WalshKieran/Integration.git
* git checkout datalake
* See very top of Integration/resources/css/style2.js
* Host it locally somehow (e.g. "web server for chrome" > choose folder > ../Integration i.e. root git folder)
* Navigate to http://localhost:8887/site.html (edited) 


# Plots

### Checklist for analysis 
- [x] Rerun all code and generate with updated pcw values
- [ ] Upload data to website 
- [ ] Move data to RNA
- [x] Update figures for BIThub
- [ ] Rerun variance parition 
- [ ] Write manual and help book for analysis - IN PROGRESS
