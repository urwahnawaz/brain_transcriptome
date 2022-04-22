## Setup
  setwd("/Volumes/Data1/PROJECTS/Urwah/integration_datasets/")

## Parameters and functions
  thresh <- 1  
  thresh.in.weak <- 0.1
  thresh.in.stringent <- 0.2
  
  apply.threshold <- function(data, t = thresh, fraction = thresh.in.weak) {
    logi <- data > t
    keep <- (rowSums(logi) / ncol(logi)) > fraction
    
    data <- data[keep,]
    return(data)
  } 
  
## Load all data
  e <- list()
  
  e$Raw$GTEx <- read.csv("FormattedData/Gtex/gtex_tpm_matrix_brainOnly.csv")
  e$Raw$PE <- read.csv("FormattedData/PsychEncode/DER-02_PEC_gene_expression_matrix_tpm_formatted.csv")
  e$Raw$BSpan <- read.csv("FormattedData/BrainSpan/brainspan_matrix_rpkm.csv")
  e$Raw$BSeq <- read.csv("FormattedData/BrainSeq/brainseq_matrix_rpkm.csv")

  e$Raw <- lapply(e$Raw, function(x) {
    rownames(x) <- x[,1]
    x <- x[,-1]
    return(x)
  })
  
## Process
  e$Scale <- lapply(e$Raw, scale) # basic z-score columnwise transformation
  e$Log <- lapply(e$Raw, function(x) { log2(x+0.5) }) # log2
  e$Log_Scale <- lapply(e$Log, scale) # log2 followed by scale
  
## Threshold, then process
  e$Thresh <- lapply(e$Raw, apply.threshold) # this is > 1 unit (TPM in GTEx, RPKM in others) in > 10% of samples
  e$Thresh_Scale <- lapply(e$Thresh, scale)
  e$Thresh_Log <- lapply(e$Thresh, function(x) { log2(x+0.5) })
  e$Thresh_Log_Scale <- lapply(e$Log, scale)
  
## Histogram of median expression
  # collect median expression values across samples for every dataset
  p <- lapply(e, function(x) {
    # lapply(x, rowMeans)
    lapply(x, function(y) {
      apply(y, 1, median)
    })
  })
  
  p <- melt(p)
  colnames(p) <- c("MedianExp", "Dataset", "Transform")
  p$Dataset <- factor(p$Dataset)
  
  # move thresholding category to a separate column
  p$Thresh <- FALSE; p$Thresh[grep("Thresh", p$Transform)] <- TRUE
  p$Transform <- gsub("Thresh_", "", p$Transform)
  p$Transform[which(p$Transform == "Thresh")] <- "Raw"
  p$Transform <- factor(p$Transform, levels = c("Raw", "Log", "Scale", "Log_Scale"))
  
  # ...and to separate dataframes
  t <- p[which(p$Thresh),]
  x <- table(t$Dataset) / 4 # this gets the number of genes in each dataset (divisor is 4 as each gene is represented by 4 transformations)
  levels(t$Dataset) <- paste0(names(x), " (n=", x, ")")
  
  p <- p[-which(p$Thresh),]
  x <- table(p$Dataset) / 4 # this gets the number of genes in each dataset (divisor is 4 as each gene is represented by 4 transformations)
  levels(p$Dataset) <- paste0(names(x), " (n=", x, ")")
  
  # plot
  pdf(file = "Expression Distributions.pdf", height = 9, width = 9)
  ggplot(p, aes(x = MedianExp, colour = Dataset)) +
    geom_density() +
    # geom_histogram() +
    facet_wrap(~Transform, scales = "free") +
    theme_bw() +
    labs(y = "Density", x = "Median Expression Across Samples", title = "No Expression Filter")
  
  ggplot(t, aes(x = MedianExp, colour = Dataset)) +
    geom_density() +
    # geom_histogram() +
    facet_wrap(~Transform, scales = "free") +
    theme_bw() +
    labs(y = "Density", x = "Median Expression Across Samples", title = ">1 RPKM in 10% of Samples")
  dev.off()
  
## Correspondence
  p <- lapply(e[c(1,4,5,8)], function(x) {
    lapply(x, function(y) {
      apply(y, 1, median)
    })
  })
  
  q <- list()
  q$Thresh <- q$Raw <- list()
  for (j in names(e$Raw)) {
    q$Raw[[j]] <- cbind(p$Raw[[j]], p$Log_Scale[[j]])
    q$Thresh[[j]] <- cbind(p$Thresh[[j]], p$Thresh_Log_Scale[[j]])
  }
  
  scatterplot <- function(x, y, lim = 100) { # x is Raw/Thresh, y is dataset
    z <- q[[x]][[y]]
    
    # filter to set limit
    rem <- (which(z[,1] > lim))
    z <- z[-rem,]
    
    # calculate correspondences
    at1 <- round(z[which.min(abs(z[,1] - 1)), 2],2)
    at10 <- round(z[which.min(abs(z[,1] - 10)), 2],2)
    at25 <- round(z[which.min(abs(z[,1] - 25)), 2],2)
    at100 <- round(z[which.min(abs(z[,1] - 100)), 2],2)
    at <- paste(at1, at10, at25, at100, sep = "/")
    
    qplot(z[,1], z[,2]) +
      theme_bw() +
      scale_x_continuous(limits = c(NA, lim)) +
      # theme(axis.title.x = element_blank(), axis.ti)
      labs(x = paste0("Untransformed RPKM/TPM", "\nx=1/10/25/100 | y=", at), 
           y = "Scaled log2+0.5", 
           title = paste(x, y, "(n=", nrow(z), ",", length(rem), "not shown)")) 
  }
  
  p <- list()
  p$AR <- scatterplot("Raw", "PE")
  p$AT <- scatterplot("Thresh", "PE")
  p$BR <- scatterplot("Raw", "BSpan")
  p$BT <- scatterplot("Thresh", "BSpan")
  p$CR <- scatterplot("Raw", "BSeq")
  p$CT <- scatterplot("Thresh", "BSeq")
  p$DR <- scatterplot("Raw", "GTEx")
  p$DT <- scatterplot("Thresh", "GTEx")

  pdf(file = "Expression Correspondence (Median).pdf", height = 16, width = 8)
  plot_grid(plotlist = p, ncol = 2)
  dev.off()
