############################################################################################################################### #
## Setup ----

## This script will perform a PCA of all four brain datasets
#setwd("/Volumes/Data1/PROJECTS/Urwah/integration_datasets/FormattedData")

## Functions
  # thresholding
    thresh <- function(x) {
      y <- x > 1
      keep <- rowSums(y) > (ncol(y) / 10)
      return(x[keep,])
    }

################################################################################################################################ #
## Load data ----    
    
## Load data
  gtex <- read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-exp.csv", row.names = 1, check.names = FALSE)
 # rownames(gtex) <- gtex$X
  #gtex <- gtex[,-1]
  m.gtex <- read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-metadata.csv")


span <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/BrainSpan-exp.csv", row.names = 1, check.names = FALSE)
rownames(span) <- span$EnsemblID
span <- span[,-1]
m.span <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/BrainSpan-metadata.csv")

seq <- read.csv("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-exp.csv")
rownames(seq) <- seq$X
seq <- seq[,-1]
m.seq <- read.csv("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-metadata.csv")

pe <- read.csv("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted/PsychEncode-exp.csv", row.names = 1, check.names = FALSE)
rownames(pe) <- pe$EnsemblID
pe <- pe[,-1]
m.pe <- read.csv("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted/PsychEncode-metadata.csv")
m.pe <- m.pe[-which(duplicated(m.pe$SampleID)),]

save(pe, seq, span, gtex, file = "/home/neuro/Documents/Brain_integrative_transcriptome/Data/Exp_Matrices_BIThub.Rda")

## Threshold
gtex <- thresh(gtex)
span <- thresh(span)
seq <- thresh(seq)
pe <- thresh(pe)

################################################################################################################################ #
## Combine data ----

## Combine into single expression matrix
  #rownames(gtex) <- sapply(strsplit(rownames(gtex), "\\."), "[", 1)  # for common formatting (this removes everything after the ".")
  common <- intersect(rownames(gtex), rownames(span))
  common <- intersect(common, rownames(seq))
  common <- intersect(common, rownames(pe))
  
  merge <- cbind(gtex[common,], 
                 span[common,], 
                 seq[common,], 
                 pe[common,])

  
  # colnames(m.pe) <- gsub("Acronym", "", colnames(m.pe))
m.pe$Regions <- "Cortex"
m.pe$StructureAcronym <- "PFC"
m.pe$DonorID <- 1:nrow(m.pe)
m.seq$StructureAcronym <- m.seq$Structure
    
    m.merge <- rbind(m.gtex[,c("SampleID", "AgeInterval", "Sex", "Period", "Regions", "StructureAcronym", "Regions")],
                   m.span[,c("SampleID", "AgeInterval", "Sex", "Period", "Regions", "StructureAcronym", "Regions")],
                   m.seq[,c( "SampleID", "AgeInterval", "Sex", "Period", "Regions", "StructureAcronym", "Regions")],
                   m.pe[,c( "SampleID", "AgeInterval", "Sex", "Period", "Regions", "StructureAcronym", "Regions")])
    m.merge$Study <- c(rep("GTEx", ncol(gtex)),
                       rep("BrainSpan", ncol(span)),
                       rep("BrainSeq", ncol(seq)),
                       rep("PsychENCODE", ncol(pe)))
    
  ## Fix age interval
    # table(m.merge$AgeInterval[which(m.merge$Study == "GTEx")])
    table(m.merge$AgeInterval)
    # first, remove spaces
    m.merge$AgeInterval <- gsub(" ", "", m.merge$AgeInterval)
    
    # add "yrs" to the GTEx samples
   # m.merge$AgeInterval[which(m.merge$Study == "GTEx")] <- paste0(m.merge$AgeInterval[which(m.merge$Study == "GTEx")], "yrs")
    
    # convert to sorted factor
    m.merge$AgeInterval <- factor(m.merge$AgeInterval, levels = levels(as.factor(m.merge$AgeInterval))[c(18, 2, 4, 5, 6, 9, 1, 15,7, 14, 3, 8, 10, 12, 13, 16, 17, 19, 20)])
    
  ## Fix period
    table(m.merge$Period)
    m.merge$Period <- gsub("p", "P", m.merge$Period)
    m.merge$Period <- factor(m.merge$Period, levels = c("Prenatal", "Postnatal")) # converts to factor and rearranges
    
################################################################################################################################ #
## PCA ----
    
## Normalise ## Non-filtered data
  # create a quantile normalised version
  library(preprocessCore)
  q.merge <- normalize.quantiles(as.matrix(merge), copy = FALSE)
  
  # log transform
  q.merge <- log2(q.merge + 0.5)
  
## PCA
  pca <- princomp(q.merge, cor = TRUE) # ~10 minutse
  
  plot.data <- as.data.frame(pca$loadings[,1:2])
  colnames(plot.data) <- c("PC1", "PC2")
  plot.data <- data.frame(plot.data, m.merge)
  
  
## Plot
  library(ggplot2)
  pdf(file = "PCA Plots.pdf", height = 6, width = 9)
  ggplot(plot.data, aes(x = PC1, y = PC2, colour = Study)) +
    geom_point()
  
  ggplot(plot.data, aes(x = PC1, y = PC2, colour = Period)) +
    geom_point()
  
  ggplot(plot.data, aes(x = PC1, y = PC2, colour = AgeInterval)) +
    geom_point()
  dev.off()
  
  
################################################################################################################################ #
## Metadata plots ----


  colnames(m.merge) 
  # "DonorID" "SampleID" "AgeInterval" "Sex" "Period" "Structure" "StructureAcronym" "Regions" "Study"
  
## Number of samples per region
  pdf(file = "Metaplots - Region.pdf", height = 4, width = 8)
  ggplot(m.merge, aes(x = Study, fill = Regions)) +
    geom_bar(colour = "black") +
    theme_bw() +
    # scale_y_continuous(expand = )
    theme(panel.border = element_blank(),
          axis.line.y = element_line()) +
    labs(y = "Sample Count")
  
  ggplot(m.merge, aes(x = Study, fill = Regions)) +
    geom_bar(colour = "black", position = "dodge") +
    theme_bw() +
    # scale_y_continuous(expand = )
    theme(panel.border = element_blank(),
          axis.line.y = element_line()) +
    labs(y = "Sample Count")
  dev.off()
  
  pdf(file = "Metaplots - Structure.pdf", height = 4, width = 8)
  ggplot(m.merge, aes(x = Study, fill = StructureAcronym)) +
    geom_bar(colour = "black") +
    theme_bw() +
    # scale_y_continuous(expand = )
    theme(panel.border = element_blank(),
          axis.line.y = element_line()) +
    labs(y = "Sample Count")
  
    ggplot(m.merge, aes(x = Study, fill = StructureAcronym)) +
    geom_bar(colour = "black", position = "dodge") +
    theme_bw() +
    # scale_y_continuous(expand = )
    theme(panel.border = element_blank(),
          axis.line.y = element_line()) +
    labs(y = "Sample Count")
  dev.off()
  
  
## Age breakdown per region
  pdf(file = "Metaplots - Age Per Region.pdf", height = 4, width = 10)
  ggplot(m.merge, aes(x = Regions, fill = AgeInterval)) +
    geom_bar() +
    theme_bw() +
    facet_wrap(~Study, nrow = 1) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Sample Count")
  
  ggplot(m.merge, aes(x = AgeInterval, fill = Regions)) +
    geom_bar() +
    theme_bw() +
    facet_wrap(~Study, nrow = 1) +
    geom_vline(xintercept = c(6.5, 11.5)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Sample Count")
  
   ggplot(m.merge, aes(x = AgeInterval, fill = Regions)) +
    geom_bar() +
    theme_bw() +
    facet_wrap(~Study, nrow = 1, scales = "free_y") +
     geom_vline(xintercept = c(6.5, 11.5)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Sample Count")
   dev.off()
  
## Total number of individuals, regions
   
   
   
                                 
