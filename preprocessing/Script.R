## Functions and libraries
source("snRNAseq-fun.R")



### INPUT DIRECTORIES 
## HCA 
hca =  file.path("/home/neuro/Documents/BrainData/single-cell/hca/Raw")
hca.outfile = file.path("/home/neuro/Documents/BrainData/single-cell/hca/Processed/")


## Velmeshev 

vel = file.path("/home/neuro/Documents/BrainData/single-cell/velmeshev/Velmeshev2019/")
vel.outfile = file.path("/home/neuro/Documents/BrainData/single-cell/velmeshev/FormattedData")

################################################################################################################################ #
## CA ----

## Read in

dat <- read.csv(file.path(hca, "matrix.csv")) # contains all brain regions
rownames(dat) <- dat$sample_name
dat <- dat[,-1]
# dat2 <- dat
dat <- t(dat)  
  
# add gene symbol
meta <- read.csv(file.path(hca, "metadata.csv"))
# dat <- dat[,-1] # remove an annotation column
# rownames(dat) <- meta$gene

# create Seurat object
obj$CA <- CreateSeuratObject(counts = dat,
                             min.cells = round(ncol(dat) / 100),
                             min.features = min.features,
                             project = "CA")


# read in metadata
# meta <- read.csv("/Volumes/Data1/PROJECTS/BrainCellularComposition/Data/Raw/human_MTG_2018-06-14_samples-columns.csv")
# obj$CA$Individual <- meta$donor
obj$CA$orig.celltype <- meta$cluster_label

## Remove cells with no class
keep <- which(!(obj$CA$orig.celltype == ""))
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
write.csv(dat, file = file.path(hca.outfile, "HCA-exp.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)

## Process metadata
  # rename cell-types
  rownames(meta) <- meta$sample_name
  meta <- meta[colnames(obj$CA),] 
  meta$MajorCelltype <- meta$cluster_label
  meta <- rename("Inh|GABA", "Inhibitory Neurons")
  meta <- rename("Exc|Glut", "Excitatory Neurons")
  meta <- rename("Astro", "Astrocytes")
  meta <- rename("Endo|VLMC|Peri", "Vasculature")
  meta <- rename("Micro", "Microglia")
  meta <- rename("Oligo", "Oligodendrocytes")
  meta <- rename("OPC", "OPCs")
  meta <- rename("Non", "Unassigned Nonneuronal")
  
  # add further donor information from the metadata of HCA release 1
  metaold <- read.csv(file.path(dir,"metadata_release1.csv"))
  m <- match(meta$external_donor_name_label, metaold$donor)
  meta$donor_age_days <- as.numeric(as.character(metaold$age_days[m]))
  meta$AgeNumeric <- meta$donor_age_days / 365
  
  # save
  
  write.csv(meta, file = file.path(hca.outfile, "HCA-metadata.csv"), quote = FALSE, row.names = FALSE)


################################################################################################################################ #
## Velmeshev ----

## Load

dat <- Read10X(vel)
obj$VL <- CreateSeuratObject(counts = dat,
                             min.cells = round(ncol(dat) / 100),
                             min.features = min.features,
                             project = "VL")
rm(dat)

## Annotate
meta <- read.table(file.path(vel, "meta.txt"), sep = "\t", header = TRUE)
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
write.csv(dat, file = file.path(vel.outfile, "Velmeshev-exp.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)

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
  write.csv(meta, file = file.path(vel.outfile, "Velmeshev-metadata.csv"), quote = FALSE, row.names = FALSE)

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


