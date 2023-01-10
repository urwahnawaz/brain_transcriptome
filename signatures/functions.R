write.gof<- function(measuredExp, estimatedComp, signatureUsed, returnPred = FALSE) {
# set to common row order
commonGenes <- rownames(measuredExp)[which(rownames(measuredExp) %in% rownames(signatureUsed))]
measuredExp <- measuredExp[commonGenes,]; signatureUsed <- signatureUsed[commonGenes,]  

# quantile normalise    
qn <- data.frame(signatureUsed, measuredExp)
qn <- as.data.frame(normalize.quantiles(as.matrix(qn), copy = FALSE))
signatureUsed <- qn[,1:ncol(signatureUsed)]
measuredExp <- qn[,-c(1:ncol(signatureUsed))]  

# predict expression (predExp) from the estimatedComp * signatureUsed
predExp <- as.data.frame(matrix(nrow = length(commonGenes), ncol = ncol(measuredExp)))
rownames(predExp) <- commonGenes    

for(j in 1:ncol(predExp)) {
    # storage
    a <- list()      
    
    # the contribution of each cell-type to predicted expression
    for(k in colnames(signatureUsed)) { a[[k]] <- estimatedComp[j,k] * signatureUsed[,k] }     
    
    # sum expression from all cell-types to a single predicted value
    predExp[,j] <- rowSums(do.call("cbind", a))
}

## Calculate statistics
stats <- as.data.frame(matrix(ncol = 5, nrow = ncol(measuredExp)))
colnames(stats) <- c("rho", "r", "mae", "rmse", "nmae")
rownames(stats) <- colnames(measuredExp)    
for(j in 1:ncol(measuredExp)) { 
    a <- measuredExp[,j]
    b <- predExp[,j]      
    stats$r[j] <- cor(log2(a+0.5), log2(b+0.5), method = "p")
    stats$rho[j] <- cor(a, b, method = "s")      
    stats$mae[j] <- mae(a, b) 
    stats$rmse[j] <- rmse(a, b)       
    stats$nmae[j] <- compute.nmae(a, b)     
}  
# return
if(returnPred) {
    res <- list()
    res$predExp <- predExp
    res$stats <- stats  
} else {
    res <- stats
}    

return(res)
}


calculate_rpkm = function(exp){
    
    ## Retrieve gene lengths from bioMart 
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
    annotations <- biomaRt::getBM(mart = ensembl, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"))
    annotations <- dplyr::transmute(annotations, ensembl_gene_id, external_gene_name, gene_length = end_position - start_position)
    
    ## Ensure gene lengths match exp matrix 
    final.genes <- annotations %>% dplyr::filter(annotations$ensembl_gene_id %in% rownames(exp))
    final.genes <- final.genes[order(match(final.genes$ensembl_gene_id, rownames(exp))),]; rownames(final.genes) <-NULL
    
    exp.rpkm = signatures[rownames(exp) %in% final.genes$ensembl_gene_id,]
    expression.rpkm <- data.frame(sapply(exp.rpkm, function(column) 10^9 * column / final.genes$gene_length / sum(column)))
    rownames(expression.rpkm) = rownames(exp.rpkm)
    return(expression.rpkm)
}