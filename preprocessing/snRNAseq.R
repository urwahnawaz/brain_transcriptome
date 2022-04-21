## In this script, mixtures are made from single-nucleus data

################################################################################################################################ #
## Setup ----

## Generic
rm(list=ls())
options(stringsAsFactors = FALSE)

## Set directory
setwd("/Volumes/Data1/PROJECTS/Urwah/integration_datasets/SingleNucleus/")

## Functions and libraries
source("/Volumes/Data1/PROJECTS/BrainCellularComposition//Scripts/Fun_Preprocessing.R")
load(file = "/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/geneInfo.rda")
load("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Preprocessed/exonicLength.rda")
library(Seurat)

## Lists
# to hold the dataset-level Seurat objects
obj <- list() # to store data from...
obj$VL <- list() # Velmeshev 2019
obj$NG <- list() # Nagy 2020
obj$CA <- list() # Hodge 2019


## Key parameters
# for seurat preprocessing
min.cells <- 0 # during the initial load, a gene is excluded if in < 0 cells 
min.features <- 200 # during the initial load, a barcode is excluded < 200 features are expressed
min.depth <- 1000 # a barcode is excluded if nCount_RNA < this value
max.depth.percentile <- 0.995 # a barcode is excluded if nCount_RNA > this percentile within the dataset
max.mito <- 5
min.celltype.n <- 0 # minimum number of members in a celltype for it to be kept. applied to anything used for creating mixtures (at this stage, Vel and HCA, but the former passes this criterion for all celltypes anyway...)

# preprocessing options
downsample <- FALSE
downsample.n <- NA; if (downsample) downsample.n <- NA
use.SCTransform <- FALSE


## Functions
## Function for downsampling the dataset to a set number of barcodes
downsample.fun <- function(x, n = downsample.n) {
  if (ncol(x) <= downsample.n) {
    print("No downsampling performed (Reason: number of cells in dataset is already less than or equal to the downsampling number)")
  } else {
    sample <- sample(colnames(x), size = n, replace = FALSE)
    x <- subset(x, cells = sample)  
  }
  return(x)
}

## General function for preprocessing sn data (normalise, filters, and scales)
get.max.depth <- function(x) {
  max.depth <- quantile(x@meta.data$nCount_RNA, probs = max.depth.percentile)
}

preprocess.fun <- function(x, run.downsample = downsample, SCTransform = use.SCTransform, max.depth = max.depth) {
  # quantify mitochondrial reads
  x[["percent.mito"]] <- PercentageFeatureSet(object = x, pattern = "^MT-")
  
  # filter to remove outlier nuclei: 
  
  x <- subset(x = x, subset = (nCount_RNA > min.depth) & (nCount_RNA < max.depth) & (percent.mito < max.mito))
  
  # downsample
  if (run.downsample) { x <- downsample.fun(x) }
  
  # normalise expression levels
  x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
  
  # find variable genes (i.e. features)
  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000)
  
  
  # further normalisation
  if (use.SCTransform) {
    x <- SCTransform(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
  }
  
  # output
  return(x)
} 

## Function for brief UMAP visualisation. x must be the output of preprocess.fun(x)
UMAP.fun <- function(x, dims = 30) {
  x <- ScaleData(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
  x <- RunPCA(x, npcs = dims)
  x <- FindNeighbors(object = x, dims = 1:dims) 
  x <- FindClusters(object = x, resolution = 1) # "We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets"
  x <- RunUMAP(object = x, dims = 1:dims)
}



## Function for library size correction
  make.cpm <- function(x) {
    for (j in 1:ncol(x)) {
      x[,j] <- x[,j] / sum(x[,j]) * 10^6
    }
    return(x)
  }

## Function to reclassify subtypes to major cell-types
rename <- function(old, new, m = meta) {
  m$MajorCelltype[grep(old, m$MajorCelltype)] <- new
  return(m)
}

################################################################################################################################ #
## CA ----

## Read in
dat <- read.csv("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/human_MTG_2018-06-14_exon-matrix.csv")

  
  
# add gene symbol
meta <- read.csv("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/human_MTG_2018-06-14_genes-rows.csv")
dat <- dat[,-1] # remove an annotation column
rownames(dat) <- meta$gene

# create Seurat object
obj$CA <- CreateSeuratObject(counts = dat,
                             min.cells = round(ncol(dat) / 100),
                             min.features = min.features,
                             project = "CA")


# read in metadata
meta <- read.csv("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/human_MTG_2018-06-14_samples-columns.csv")
obj$CA$Individual <- meta$donor
obj$CA$orig.celltype <- meta$cluster



## Remove cells with no class
keep <- which(!(obj$CA$orig.celltype == "no class"))
obj$CA <- subset(obj$CA, cells = keep)

## Preprocess
max.depth <- get.max.depth(obj$CA)
obj$CA <- preprocess.fun(obj$CA, max.depth = max.depth)


## Save CPM
dat <- as.data.frame(obj$CA@assays$RNA@counts)
dat <- make.cpm(dat)
dat <- cbind(rownames(dat), dat)
colnames(dat)[1] <- "Symbol"
rownames(dat) <- 1:nrow(dat)
write.csv(dat, file = "FormattedData/HCA-exp.csv", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Process metadata
  # rename cell-types
  meta$MajorCelltype <- meta$cluster
  meta <- rename("Inh|GABA", "Inhibitory Neurons")
  meta <- rename("Exc|Glut", "Excitatory Neurons")
  meta <- rename("Astro", "Astrocytes")
  meta <- rename("Endo", "Endothelia")
  meta <- rename("Micro", "Microglia")
  meta <- rename("Oligo", "Oligodendrocytes")
  meta <- rename("OPC", "OPCs")
  meta <- rename("Non", "Unassigned Nonneuronal")
  
  rownames(meta) <- meta$sample_name
  meta <- meta[colnames(obj$CA),] 
  write.csv(meta, file = "FormattedData/HCA-metadata.csv", quote = FALSE, row.names = FALSE)


################################################################################################################################ #
## Velmeshev ----

## Load
dat <- Read10X("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/Velmeshev2019/")
obj$VL <- CreateSeuratObject(counts = dat,
                             min.cells = round(ncol(dat) / 100),
                             min.features = min.features,
                             project = "VL")
rm(dat)

## Annotate
meta <- read.table("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/Velmeshev2019/meta.txt", sep = "\t", header = TRUE)
m <- match(colnames(obj$VL), meta$cell)

obj$VL$orig.celltype <- meta$cluster[m]
obj$VL$Region <- meta$region[m]
obj$VL$Disorder <- meta$diagnosis[m]
obj$VL$orig.ident <- paste0("Vel_", substr(meta$cell[m], start = 1, stop = 16))
obj$VL$Individual <- meta$individual[m]

## Preprocess and normalise
max.depth <- get.max.depth(obj$VL)
obj$VL <- preprocess.fun(obj$VL, max.depth = max.depth)

## CPM
dat <- as.data.frame(obj$VL@assays$RNA@counts)
dat <- make.cpm(dat)
dat <- cbind(rownames(dat), dat)
colnames(dat)[1] <- "Symbol"
rownames(dat) <- 1:nrow(dat)
write.csv(dat, file = "FormattedData/Velmeshev-exp.csv", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Save metadata
  # rename celltypes 
  meta$MajorCelltype <- meta$cluster
  meta <- rename("IN", "Inhibitory Neurons")
  meta <- rename("^L", "Excitatory Neurons")
  meta <- rename("AST", "Astrocytes")
  meta <- rename("Endo", "Endothelia")
  meta <- rename("Micro", "Microglia")
  meta <- rename("Oligo", "Oligodendrocytes")
  meta <- rename("OPC", "OPCs")
  meta <- rename("NRGN", "NRGN Neurons")
  meta <- rename("mat", "Maturing Neurons")

  # save
  m <- match(colnames(obj$VL), meta$cell)
  meta <-  meta[m,]
  write.csv(meta, file = "FormattedData/Velmeshev-metadata.csv", quote = FALSE, row.names = FALSE)

################################################################################################################################ #
## Nagy ----

## Load 
dat <- Read10X("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/Nagy2020/")
obj$NG <- CreateSeuratObject(counts = dat,
                             min.cells = round(ncol(dat) / 100),
                             min.features = min.features,
                             project = "Nagy")
rm(dat)

## Augment the metadata
meta <- strsplit(colnames(obj$NG), "\\.")
meta <- lapply(meta, function(x) {
  x <- c(x[1], strsplit(x[2], "_")[[1]])
  return(x)
})

meta <- do.call("rbind", meta)
meta <- as.data.frame(meta)
colnames(meta) <- c("Celltype", 
                    "Individual", # inferred from there being 17 CTL and 17 MDD patients in the study, and there are 17 levels that are always Control and 17 always Suicide in $Disorder 
                    "Disorder", # note that suicide is equivalent to MDD in this study
                    "Batch", # likely batch given its correspondence across individuals, and that it has six levels which is the number of reported batches. not used by us, so of no importance
                    "Barcode")

obj$NG$orig.celltype <- meta$Celltype
obj$NG$orig.ident <- paste0("Nagy_", meta$Barcode)
obj$NG$Individual <- meta$Individual
obj$NG$Disorder <- meta$Disorder

## Minor cell filtering: remove cells labelled as "Mix", or those from MDD individuals
# keep1 <- which(obj$NG$Disorder == "Control")
keep <- grep("Mix_", obj$NG$orig.celltype, invert = TRUE)
# keep <- intersect(keep1, keep2)
obj$NG <- subset(obj$NG, cells = keep) 

## Preprocess and normalise
max.depth <- get.max.depth(obj$NG)
obj$NG <- preprocess.fun(obj$NG, max.depth = max.depth)

## CPM
dat <- as.data.frame(obj$NG@assays$RNA@counts)
dat <- make.cpm(dat)
dat <- cbind(rownames(dat), dat)
colnames(dat)[1] <- "Symbol"
rownames(dat) <- 1:nrow(dat)
write.csv(dat, file = "FormattedData/Nagy-exp.csv", row.names = TRUE, col.names = TRUE, quote = FALSE)


