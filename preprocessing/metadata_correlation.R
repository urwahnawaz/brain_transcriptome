# Corr plots 
# Loading all libraries
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

md.df = read.csv("/home/neuro/Documents/Coc/datasets/FormattedData/FormattedData/Gtex/GTEx-metadata.csv", row.names = 1) %>%
  as.data.frame()

M = data.matrix(md.df)

M = cor(M)

corrplot(M)
M


M = cor(data.matrix(final), use = "complete.obs")
pdf("../../../Brain_integrative_transcriptome/Results/Correlations/GTEx-metadata.pdf", height = 20, width = 25)
corrplot(M)
dev.off()
<<<<<<< HEAD
=======




>>>>>>> main
