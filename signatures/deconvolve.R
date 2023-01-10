## dtangle 
source("libraries.R")


bspan.exp = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-exp.csv", row.names = 2, check.names = FALSE)[,-1]
bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv", row.names = 1, header= TRUE)


prenatal = bspan.md %>% dplyr::filter(Period == "Prenatal")
prenatal.exp = bspan.exp %>% dplyr::select(prenatal$SampleID)

signatures %<>% dplyr::select(-c("Micro.Neonatal", "Vas.Fetal", "Vas.Neonatal"))

common_genes = intersect(rownames(prenatal.exp), rownames(expression.rpkm))
prenatal.exp <- prenatal.exp[pmatch(common_genes, rownames(prenatal.exp)), ]
expression.rpkm <- expression.rpkm[pmatch(common_genes, rownames(expression.rpkm)), ]

y <- cbind(expression.rpkm, prenatal.exp) %>% dplyr::select(-c("Micro.Neonatal", "Vas.Fetal", "Vas.Neonatal"))
y <- normalizeBetweenArrays(y)
y <- t(y)



annot = colnames(signatures) %>% 
    as.data.frame() %>% 
    set_colnames("Sample") %>% 
    mutate(Cell_type = gsub("\\..*", "", Sample))

all_ct = unique(annot$Cell_type)

all_ct

head(annot)

annot$Cell_type

ps = lapply(1:length(all_ct), function(i) {
    which(annot$Cell_type == all_ct[i])
})

names(ps) = all_ct



marker_list = find_markers(y,pure_samples=ps,data_type="rna-seq",marker_method='ratio')

q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(ps)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_markers


marks = marker_list$L
dc <- dtangle(y, pure_samples=ps, n_markers=n_markers, data_type = 'rna-seq', markers = marks)
final_est <- dc$estimates[(dim(signatures)[2]+1):dim(y)[1],]
colnames(final_est) <-  all_ct

head(final_est)


plot_data <- melt(final_est)
colnames(plot_data) <- c("Sample", "Cell Type", "Proportion") 

plot_data$Proportion <- as.numeric(plot_data$Proportion)

ggplot(plot_data, aes(x = `Cell Type`, y=Proportion))+geom_violin(aes(fill = `Cell Type`)) + geom_jitter(height = 0, width = 0.1)