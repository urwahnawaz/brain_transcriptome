## dtangle 
source("libraries.R")
source("functions.R")

load("../../Results/signatures/pfc_signatures.Rda")
load("../../DeconRNAShiny/sigsBrain.rda")

bspan.exp = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-exp.csv", row.names = 2, check.names = FALSE)[,-1]
bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv", row.names = 1, header= TRUE)


bspan.md$AgeInterval
prenatal = bspan.md %>% dplyr::filter(Period == "Prenatal" & Regions == "Cortex")

prenatal.exp = bspan.exp %>% dplyr::select(prenatal$SampleID)

pfc_signatures$rpkm_all_neuro

sigs_all_prenatal = pfc_signatures$rpkm_all_neuro %>% 
    dplyr::select(contains("Fetal"))

sigs_all_adult = pfc_signatures$rpkm_all_neuro %>% 
    dplyr::select(contains("Adult"))

res = list()


sigs_all_prenatal

res$prenatal =run_dtg(bspan.exp, sigs_all_prenatal) %>% 
    as.data.frame() %>% 
    mutate(Sig = "Prenatal")

res$adult = run_dtg(bspan.exp, sigs_all_adult) %>% 
    as.data.frame() %>% 
    mutate(Sig = "Adult")

res$prenatal_and_adult = run_dtg(bspan.exp, pfc_signatures$rpkm_all_neuro_fetal_adult) %>% 
    as.data.frame() %>% 
    mutate(Sig = "Prenatal and Adult")




## subset prenatal samples 
prenatal = bspan.md %>% dplyr::filter(Period == "Prenatal")

prenatal.exp = bspan.exp %>% dplyr::select(prenatal$SampleID)
res_prenatal = res$prenatal %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    dplyr::filter(Sample %in% colnames(prenatal.exp)) %>% 
    dplyr::select(-contains("Sig"))

res_prenatal_adult = res$prenatal_and_adult %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    dplyr::filter(Sample %in% colnames(prenatal.exp)) %>% 
    dplyr::select(-contains("Sig"))
    


## GoF 

gof_res = write.gof(prenatal.exp, res_prenatal, 
                    signatureUsed = sigs_all_prenatal)


gof_both_sigs = write.gof(prenatal.exp, res_prenatal_adult, 
                    signatureUsed = pfc_signatures$rpkm_all_neuro_fetal_adult)
head(gof_res)

gof_res %>% 
    mutate(col = "Prenatal Signatures") %>%
    ggplot(aes(col, r)) + geom_violin()


gof_both_sigs %>% 
    mutate(col = "Prenatal and adult") %>%
    ggplot(aes(col, r)) + geom_violin()

## How does multibrain compete with this?
load("../../DeconRNAShiny/sigsBrain.rda")

calculate_rpkm(sigsBrain$MB) %>% 
    .[which(apply(., 1, max) > exp_thresh),] 
    

######## 
exp = bspan.exp 
sig = sigs_all_prenatal


## Get Length 

gene_list = getLength(rownames(sigsBrain$MB))

Msigs = sigsBrain$MB %>% 
    .[rownames(.) %in% gene_list$ensembl_gene_id,] %>%
    rpkm(., gene_list$gene_length) %>% 
    .[which(apply(., 1, max) > exp_thresh),]


z


   
res$mBrain_raw = run_dtg(bspan.exp,sigsBrain$MB)


res_mBrain_raw = res$mBrain_raw  %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    dplyr::filter(Sample %in% colnames(prenatal.exp))

gof_mSig_raw = write.gof(prenatal.exp, res_mBrain_raw, 
                          signatureUsed = sigsBrain$MB)


gof_mSig_raw %>% 
    mutate(col = "MB") %>%
    ggplot(aes(col, r)) + geom_violin()


sigsBrain$MB



### run deconvolution in these steps 


stages = c("Fetal", "Neonatal", "Infancy","Adolescence", "Adult")

dtg_res = list()
for (stage in stages){
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    
    dev_res = run_dtg(bspan.exp,sigs_dev) %>%
        as.data.frame() %>% 
        mutate(Sig = stage)
    
    dtg_res[[paste0(stage, "_PFC_sig")]] = dev_res
}


## 2) run dtangle for MB

dtg_res$MB = run_dtg(bspan.exp,sigsBrain$MB) %>% 
    as.data.frame() %>%
    mutate(Sig = MB)



## gof prenatal 


stages = c("Fetal", "Neonatal", "Infancy","Adolescence", "Adult")

dtg_res = list()
for (stage in stages){
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    
    dev_res = run_dtg(bspan.exp,sigs_dev) %>%
        as.data.frame() %>% 
        mutate(Sig = stage)
    
    dtg_res[[paste0(stage)]] = dev_res
    
}



dtg_res$MB = run_dtg(bspan.exp,sigsBrain$MB) %>% 
    as.data.frame() %>%
    mutate(Sig = MB)



## GoF for Lister Signature 
## GoF for Lister Signature 
stages = c("Fetal", "Adult")
gof_res = list()

for (stage in stages){
    if (stage == "Fetal") {
        stage.md = bspan.md %>% dplyr::filter(Period == "Prenatal")
        
    } else {
        stage.md = bspan.md %>% dplyr::filter(Period == "Postnatal")
    }
    
    stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
    dtg.stage =  dtg_res[[paste0(stage)]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    gof = write.gof(stage.exp, dtg.stage, sigs_dev) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_PFC_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    gof_res[[paste0(stage, "_gof_res")]] = gof
    
    ## MultiBrain 
    
    dtg.MB = dtg_res[["MB"]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    
    gof.MB = write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_MB_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    
    gof_res[[paste0(stage, "_gof_MB")]] = gof.MB
    
}

dtg_res$MB = run_dtg(bspan.exp,sigsBrain$MB)

dtg_res$MB %<>% as.data.frame() %>%
    mutate(Sig = MB)

dtg_res[["MB"]] = run_dtg(bspan.exp,sigsBrain$MB) %>% 
    as.data.frame() %>%
    mutate(Sig = MB)




?left_join


head(gof_res$Fetal_gof_res)
gof_res$Adult_gof_res %>% 
    ggplot(aes(Regions,r, fill = Regions)) + geom_boxplot()

gof_res$Fetal_gof_res %>% 
    ggplot(aes(Stage,r)) + geom_violin()

head(gof_res$Fetal_gof_res)
## correlation 
gof_res$Fetal_gof_res %>% 
    ggplot(aes(x=r,y=mRIN)) + geom_point()


### - deconRNAseq
library(DeconRNASeq)




## Rerun this in a function 



stages = c("Fetal", "Neonatal", "Infancy","Adolescence", "Adult")

decon_res = list()
for (stage in stages){
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    
    dev_res = run_decon(bspan.exp,sigs_dev) %>%
        as.data.frame() %>% 
        mutate(Sig = stage)
    
    decon_res[[paste0(stage, "_PFC_sig")]] = dev_res
}


decon_res$MB = run_decon(bspan.exp, sigsBrain$MB) %>% 
    as.data.frame() %>% 
    mutate(Sig = "MB")



decon_res$Fetal_PFC_sig

## Gof 

stages = c("Fetal", "Adult")
gof_res_decon = list()

for (stage in stages){
    if (stage == "Fetal") {
        stage.md = bspan.md %>% dplyr::filter(Period == "Prenatal")
        
    } else {
        stage.md = bspan.md %>% dplyr::filter(Period == "Postnatal")
    }
    
    stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
    decon.stage =  decon_res[[paste0(stage,"_PFC_sig")]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    gof = write.gof(stage.exp, decon.stage, sigs_dev) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_PFC_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    gof_res_decon[[paste0(stage, "_gof_res")]] = gof
    ## MB
    decon.MB = decon_res$MB %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    
    gof.MB = write.gof(stage.exp, decon.MB, sigsBrain$MB) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_MB_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    
    gof_res_decon[[paste0(stage, "_gof_MB")]] = gof.MB
    
    
}



stage.md = stage.md = bspan.md %>% dplyr::filter(Stage == "s11" | Stage == "s12" | Stage == "s13")
stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
decon.stage =  decon_res[[paste0("Adult","_PFC_sig")]] %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample") %>% 
    dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
    column_to_rownames("Sample") %>% 
    dplyr::select(-contains(".Sig"))
sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
    dplyr::select(contains("Adult"))


gof.ad = write.gof(stage.exp, decon.stage, sigs_dev) %>% 
    as.data.frame() %>% 
    mutate(sigs = paste0("Adult", "PFC_sig")) %>% 
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)


gof_res_decon$Adult_gof_res%>% 
    rbind(., gof_res_decon$Adult_gof_MB, gof, gof.ad) %>%
    ggplot(aes(sigs, r, color = Stage)) + geom_violin() +
    geom_jitter() + theme_bw() +  ggtitle("Goodness of fit for prenatal samples") +
    ylab("Goodness of fit (r)") + xlab("Signatures")





### create additional signature: prenatal and postnatal (all postnatal samples summed)
stage.exp 

bspan.md %>%
    dplyr::filter(Stage == "s8")

gof_res_decon$Fetal_gof_res

gof_res_decon$Fetal_gof_res %>% 
    rbind(., gof_res_decon$Fetal_gof_MB) %>%
    ggplot(aes(sigs, r)) + geom_violin() +
    geom_jitter() + theme_bw() + 
    theme(legend.position = "none") + ggtitle("Goodness of fit for prenatal samples") +
    ylab("Goodness of fit (r)") + xlab("Signatures")

gof_res_decon$Adult_gof_res%>% 
    rbind(., gof_res_decon$Adult_gof_MB) %>%
    ggplot(aes(sigs, r, color = Stage)) + geom_violin() +
    geom_jitter() + theme_bw() +  ggtitle("Goodness of fit for prenatal samples") +
    ylab("Goodness of fit (r)") + xlab("Signatures")

decon_res$Neonatal_PFC_sig
gof_res_decon
stage= "Fetal"

decon_res[[paste0(stage, "_PFC_sig")]]
decon_res



### DeconRNAseq -recommended by Gavin 
stages = c("Fetal", "Neonatal", "Infancy","Adolescence", "Adult")

decon_res = list()
for (stage in stages){
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    
    dev_res = run_decon(bspan.exp,sigs_dev) %>%
        as.data.frame() %>% 
        mutate(Sig = stage)
    
    decon_res[[paste0(stage, "_PFC_sig")]] = dev_res
}


decon_res$MB = run_decon(bspan.exp, sigsBrain$MB) %>% 
    as.data.frame() %>% 
    mutate(Sig = "MB")


max_est = apply(decon_res$Fetal[,-6],1, max)
max_est

hist(max_est)

#decon_res$Fetal_PFC_sigmax_est_fetal = apply(dtg_res$Fetal[,-6],1, max)


stages.df = list()
stages.df[["Fetal"]] = c(paste(c(4:40), "pcw", sep = "_"))
stages.df[["Neonatal"]] =c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Childhood"]] = paste(c(1:10), "yrs", sep = "_")
stages.df[["Adolescence"]] = paste(c(11:19), "yrs", sep = "_")
stages.df[["Adult"]] = paste(c(20:40), "yrs", sep = "_")


stages.df
add_feature = function(feature_column, features){
    as.vector(sapply(feature_column, function(x){
        names(features)[sapply(features, function(f) x %in% f)]})) 
}

bspan.md = bspan.md %>% 
    mutate(Dev.stage = add_feature(.$Age, stages.df))

table(bspan.md$Dev.stage)

## Goodness of fit
## Gof 

### for each deconvolution stage

## want to see well each stage does each matching stage 

### add an additional column in bspan md matching with stages 

### then want to see how each postnatal stage deconvolutes all postnatal samples well 

## Infancy in this one
stages = c("Fetal", "Neonatal","Adolescence", "Adult")


gof_res_decon = list() 

for (stage in stages){
    stage.md = bspan.md %>% dplyr::filter(Dev.stage == stage)
    stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
    decon.stage =  decon_res[[paste0(stage,"_PFC_sig")]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    gof = write.gof(stage.exp, decon.stage, sigs_dev) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_PFC_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    gof_res_decon[[paste0(stage, "_gof_res")]] = gof
    ## MB
    decon.MB = decon_res$MB %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    
    gof.MB = write.gof(stage.exp, decon.MB, sigsBrain$MB) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_MB_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    
    gof_res_decon[[paste0(stage, "_gof_MB")]] = gof.MB
}



 gof_res_decon %>% do.call(rbind, .) %>% 
     ggplot(aes(sigs, r)) + geom_violin() +facet_wrap(~Dev.stage) + 
     theme(axis.text.x = element_text(angle = 90))

gof_adult_prenatal = list()


stages = c("Fetal", "Adolescence")
stages



for (stage in stages){
    if (stage == "Fetal") {
        stage.md = bspan.md %>% dplyr::filter(Period == "Prenatal")
        
    } else {
        stage.md = bspan.md %>% dplyr::filter(Period == "Postnatal")
    }
    
    stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
    decon.stage =  decon_res[[paste0(stage,"_PFC_sig")]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    gof = write.gof(stage.exp, decon.stage, sigs_dev) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_PFC_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    gof_adult_prenatal[[paste0(stage, "_gof_res")]] = gof
    ## MB
    decon.MB = decon_res$MB %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    
    gof.MB = write.gof(stage.exp, decon.MB, sigsBrain$MB) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_MB_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    
    gof_adult_prenatal[[paste0(stage, "_gof_MB")]] = gof.MB
    
    
}



gof_adult_prenatal %>% do.call(rbind, .) %>% 
    ggplot(aes(sigs, r)) + geom_violin() +facet_wrap(~Period) + 
    theme(axis.text.x = element_text(angle = 90))


gof_adult_prenatal$Childhood_gof_res
## create two signatures 

## neonatal 
table(bspan.md$Age)


######### Deconvolve here 

source("libraries.R")
source("functions.R")

load("../../Results/signatures/pfc_signatures.Rda")
load("../../DeconRNAShiny/sigsBrain.rda")

bspan.exp = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-exp.csv", row.names = 2, check.names = FALSE)[,-1]
bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv", row.names = 1, header= TRUE)


stages.df = list()
stages.df[["Fetal"]] = c(paste(c(4:40), "pcw", sep = "_"))
#stages.df[["Neonatal"]] =c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Infancy"]] =c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Childhood"]] = paste(c(1:10), "yrs", sep = "_")
stages.df[["Adolescence"]] = paste(c(11:19), "yrs", sep = "_")
stages.df[["Adult"]] = paste(c(20:40), "yrs", sep = "_")


add_feature = function(feature_column, features){
    as.vector(sapply(feature_column, function(x){
        names(features)[sapply(features, function(f) x %in% f)]})) 
}

bspan.md = bspan.md %>% 
    mutate(Dev.stage = add_feature(.$Age, stages.df))

table(bspan.md$Dev.stage)



stages = c("Fetal", "Neonatal","Infancy",  "Childhood","Adolescence", "Adult")

dtg_res = list()
for (stage in stages){
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    
    dev_res = run_dtg(bspan.exp,sigs_dev) %>%
        as.data.frame() %>% 
        mutate(Sig = stage)
    
    dtg_res[[paste0(stage)]] = dev_res
}




## 2) run dtangle for MB

dtg_res$MB = run_dtg(bspan.exp,sigsBrain$MB) %>%
    as.data.frame() %>%
    mutate(Sig = "MB")

stages = c("Fetal", "Infancy", "Childhood","Adolescence", "Adult")
gof_res = list()

for (stage in stages){
    stage.md = bspan.md %>% dplyr::filter(Dev.stage == stage)
    stage.exp = bspan.exp %>% dplyr::select(stage.md$SampleID)
    dtg.stage =  dtg_res[[paste0(stage)]] %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    sigs_dev = pfc_signatures$rpkm_all_neuro %>% 
        dplyr::select(contains(stage))
    gof = write.gof(stage.exp, dtg.stage, sigs_dev) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_PFC_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    gof_res[[paste0(stage, "_gof_res")]] = gof
    
    ## MultiBrain 
    
    dtg.MB = dtg_res$MB %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample") %>% 
        dplyr::filter(Sample %in% colnames(stage.exp)) %>% 
        column_to_rownames("Sample") %>% 
        dplyr::select(-contains(".Sig"))
    
    gof.MB = write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>% 
        as.data.frame() %>% 
        mutate(sigs = paste0(stage, "_MB_sig")) %>% 
        rownames_to_column("SampleID") %>%
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>% 
        mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
        left_join(stage.md)
    
    gof_res[[paste0(stage, "_gof_MB")]] = gof.MB
    
}




gof_res %>% 
    do.call(rbind, .) %>% 
    ggplot(aes(sigs,r, fill = sigs)) + geom_violin() + facet_wrap(~Dev.stage) + 
    theme(axis.text.x = element_blank())




