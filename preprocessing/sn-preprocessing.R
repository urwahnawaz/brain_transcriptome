## single-cell preprocessing 


## Speir et al 
path = file.path("/home/neuro/Documents/BrainData/single-cell/speir")
md = read.table(file.path(path, "meta.tsv"), sep = "\t", header=TRUE)
md
table(md$Cluster)


