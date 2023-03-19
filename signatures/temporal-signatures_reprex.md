# Creating temporal deconvolution signatures

``` r
source("libraries.R")
#;-) Loading required package: limma
#;-) 
#;-) Attaching package: 'dplyr'
#;-) The following objects are masked from 'package:stats':
#;-) 
#;-)     filter, lag
#;-) The following objects are masked from 'package:base':
#;-) 
#;-)     intersect, setdiff, setequal, union
#;-) 
#;-) Attaching package: 'data.table'
#;-) The following objects are masked from 'package:dplyr':
#;-) 
#;-)     between, first, last
#;-) ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
#;-) ✔ ggplot2 3.4.0     ✔ purrr   1.0.1
#;-) ✔ tibble  3.1.8     ✔ stringr 1.5.0
#;-) ✔ tidyr   1.3.0     ✔ forcats 0.5.2
#;-) ✔ readr   2.1.3     
#;-) ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#;-) ✖ data.table::between() masks dplyr::between()
#;-) ✖ dplyr::filter()       masks stats::filter()
#;-) ✖ data.table::first()   masks dplyr::first()
#;-) ✖ dplyr::lag()          masks stats::lag()
#;-) ✖ data.table::last()    masks dplyr::last()
#;-) ✖ purrr::transpose()    masks data.table::transpose()
#;-) 
#;-) Attaching package: 'magrittr'
#;-) 
#;-) 
#;-) The following object is masked from 'package:purrr':
#;-) 
#;-)     set_names
#;-) 
#;-) 
#;-) The following object is masked from 'package:tidyr':
#;-) 
#;-)     extract
#;-) 
#;-) 
#;-) ------------------------------------------------------------------------------
#;-) 
#;-) You have loaded plyr after dplyr - this is likely to cause problems.
#;-) If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
#;-) library(plyr); library(dplyr)
#;-) 
#;-) ------------------------------------------------------------------------------
#;-) 
#;-) 
#;-) Attaching package: 'plyr'
#;-) 
#;-) 
#;-) The following object is masked from 'package:purrr':
#;-) 
#;-)     compact
#;-) 
#;-) 
#;-) The following objects are masked from 'package:dplyr':
#;-) 
#;-)     arrange, count, desc, failwith, id, mutate, rename, summarise,
#;-)     summarize
#;-) 
#;-) 
#;-) Loading required package: BiocGenerics
#;-) 
#;-) 
#;-) Attaching package: 'BiocGenerics'
#;-) 
#;-) 
#;-) The following objects are masked from 'package:dplyr':
#;-) 
#;-)     combine, intersect, setdiff, union
#;-) 
#;-) 
#;-) The following object is masked from 'package:limma':
#;-) 
#;-)     plotMA
#;-) 
#;-) 
#;-) The following objects are masked from 'package:stats':
#;-) 
#;-)     IQR, mad, sd, var, xtabs
#;-) 
#;-) 
#;-) The following objects are masked from 'package:base':
#;-) 
#;-)     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#;-)     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#;-)     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#;-)     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#;-)     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#;-)     union, unique, unsplit, which.max, which.min
#;-) 
#;-) 
#;-) Loading required package: S4Vectors
#;-) 
#;-) Loading required package: stats4
#;-) 
#;-) 
#;-) Attaching package: 'S4Vectors'
#;-) 
#;-) 
#;-) The following object is masked from 'package:plyr':
#;-) 
#;-)     rename
#;-) 
#;-) 
#;-) The following object is masked from 'package:tidyr':
#;-) 
#;-)     expand
#;-) 
#;-) 
#;-) The following objects are masked from 'package:data.table':
#;-) 
#;-)     first, second
#;-) 
#;-) 
#;-) The following objects are masked from 'package:dplyr':
#;-) 
#;-)     first, rename
#;-) 
#;-) 
#;-) The following objects are masked from 'package:base':
#;-) 
#;-)     expand.grid, I, unname
#;-) 
#;-) 
#;-) Loading required package: IRanges
#;-) 
#;-) 
#;-) Attaching package: 'IRanges'
#;-) 
#;-) 
#;-) The following object is masked from 'package:plyr':
#;-) 
#;-)     desc
#;-) 
#;-) 
#;-) The following object is masked from 'package:purrr':
#;-) 
#;-)     reduce
#;-) 
#;-) 
#;-) The following object is masked from 'package:data.table':
#;-) 
#;-)     shift
#;-) 
#;-) 
#;-) The following objects are masked from 'package:dplyr':
#;-) 
#;-)     collapse, desc, slice
#;-) 
#;-) 
#;-) Loading required package: GenomeInfoDb
#;-) 
#;-) Loading required package: GenomicRanges
#;-) 
#;-) 
#;-) Attaching package: 'GenomicRanges'
#;-) 
#;-) 
#;-) The following object is masked from 'package:magrittr':
#;-) 
#;-)     subtract
#;-) 
#;-) 
#;-) Loading required package: AnnotationDbi
#;-) 
#;-) Loading required package: Biobase
#;-) 
#;-) Welcome to Bioconductor
#;-) 
#;-)     Vignettes contain introductory material; view with
#;-)     'browseVignettes()'. To cite Bioconductor, see
#;-)     'citation("Biobase")', and for packages 'citation("pkgname")'.
#;-) 
#;-) 
#;-) 
#;-) Attaching package: 'AnnotationDbi'
#;-) 
#;-) 
#;-) The following object is masked from 'package:dplyr':
#;-) 
#;-)     select
#;-) 
#;-) 
#;-) Loading required package: MatrixGenerics
#;-) 
#;-) Loading required package: matrixStats
#;-) 
#;-) 
#;-) Attaching package: 'matrixStats'
#;-) 
#;-) 
#;-) The following objects are masked from 'package:Biobase':
#;-) 
#;-)     anyMissing, rowMedians
#;-) 
#;-) 
#;-) The following object is masked from 'package:plyr':
#;-) 
#;-)     count
#;-) 
#;-) 
#;-) The following object is masked from 'package:dplyr':
#;-) 
#;-)     count
#;-) 
#;-) 
#;-) 
#;-) Attaching package: 'MatrixGenerics'
#;-) 
#;-) 
#;-) The following objects are masked from 'package:matrixStats':
#;-) 
#;-)     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#;-)     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#;-)     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#;-)     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#;-)     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#;-)     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#;-)     colWeightedMeans, colWeightedMedians, colWeightedSds,
#;-)     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#;-)     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#;-)     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#;-)     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#;-)     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#;-)     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#;-)     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#;-)     rowWeightedSds, rowWeightedVars
#;-) 
#;-) 
#;-) The following object is masked from 'package:Biobase':
#;-) 
#;-)     rowMedians
#;-) 
#;-) 
#;-) 
#;-) Attaching package: 'ModelMetrics'
#;-) 
#;-) 
#;-) The following object is masked from 'package:edgeR':
#;-) 
#;-)     gini
#;-) 
#;-) 
#;-) The following object is masked from 'package:base':
#;-) 
#;-)     kappa
#;-) 
#;-) 
#;-) Loading required package: limSolve
#;-) 
#;-) 
#;-) Attaching package: 'limSolve'
#;-) 
#;-) 
#;-) The following object is masked from 'package:ggplot2':
#;-) 
#;-)     resolution
#;-) 
#;-) 
#;-) Loading required package: pcaMethods
#;-) 
#;-) 
#;-) Attaching package: 'pcaMethods'
#;-) 
#;-) 
#;-) The following object is masked from 'package:stats':
#;-) 
#;-)     loadings
#;-) 
#;-) 
#;-) Loading required package: grid
#;-) 
#;-) 
#;-) Attaching package: 'DeconRNASeq'
#;-) 
#;-) 
#;-) The following object is masked from 'package:ModelMetrics':
#;-) 
#;-)     rmse
source("functions.R")
```

``` r

## signatures
load("../../Results/signatures/pfc_signatures.Rda")
load("../../DeconRNAShiny/sigsBrain.rda")


## Brain exp
bspan.exp <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-exp.csv", row.names = 2, check.names = FALSE)[, -1]
bspan.md <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv", row.names = 1, header = TRUE)
```

``` r

## 1) subset signature for each stage and then run dtangle
stages <- c("Fetal", "Neonatal", "Infancy", "Adolescence", "Adult")

dtg_res <- list()
for (stage in stages) {
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))

  dev_res <- run_dtg(bspan.exp, sigs_dev) %>%
    as.data.frame() %>%
    mutate(Sig = stage)

  dtg_res[[paste0(stage)]] <- dev_res
}
#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf



## 2) run dtangle for MB

dtg_res$MB <- run_dtg(bspan.exp, sigsBrain$MB) %>%
  as.data.frame() %>%
  mutate(Sig = "MB")
#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf
```

``` r

## estimation vs good of fit - TO DO
### Implies deconvolution isn't good
## x axis (max of the different cell-types)
stage.md <- bspan.md %>% dplyr::filter(Period == "Prenatal")

max_est_fetal <- apply(dtg_res$Fetal[, -6], 1, max)

max_est_fetal %<>% melt() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% stage.md$SampleID) %>%
  left_join(gof_res$Fetal_gof_res)
#;-) Warning in melt(.): The melt generic in data.table has been passed a numeric
#;-) and will attempt to redirect to the relevant reshape2 method; please note that
#;-) reshape2 is deprecated, and this redirection is now deprecated as well. To
#;-) continue using melt methods from reshape2 while both libraries are attached,
#;-) e.g. melt.list, you can prepend the namespace like reshape2::melt(.). In the
#;-) next version, this warning will become an error.
#;-) Error in is.data.frame(y): object 'gof_res' not found

m <- match(gof_res$Fetal_gof_res$SampleID, names(max_est_fetal))
#;-) Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'match': object 'gof_res' not found

cor(max_est_fetal[m], gof_res$Fetal_gof_res$r, use = "pairwise.complete.obs") ## weak correlation
#;-) Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'cor': object 'm' not found

## something wrong with the deconvolution
### deconRNAseq - try
```

## Goodness of fit for just fetal and adult samples (matched)

``` r
## GoF for Lister Signature
stages <- c("Fetal", "Adult")
gof_res <- list()

for (stage in stages) {
  if (stage == "Fetal") {
    stage.md <- bspan.md %>% dplyr::filter(Period == "Prenatal")
  } else {
    stage.md <- bspan.md %>% dplyr::filter(Period == "Postnatal")
  }

  stage.exp <- bspan.exp %>% dplyr::select(stage.md$SampleID)
  dtg.stage <- dtg_res[[paste0(stage)]] %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))
  gof <- write.gof(stage.exp, dtg.stage, sigs_dev) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_PFC_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)
  gof_res[[paste0(stage, "_gof_res")]] <- gof

  ## MultiBrain

  dtg.MB <- dtg_res$MB %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))

  gof.MB <- write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_MB_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)

  gof_res[[paste0(stage, "_gof_MB")]] <- gof.MB
}
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
```

``` r

gof_res$Fetal_gof_res %>%
  rbind(., gof_res$Fetal_gof_MB) %>%
  ggplot(aes(sigs, r, fill = sigs)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).
#;-) Warning: Removed 10 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/wgeDypC.png)<!-- -->

``` r
gof_res$Fetal_gof_res %>%
  ggplot(aes(Regions, r, fill = Regions)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples, correlation with regions") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 5 rows containing non-finite values (`stat_ydensity()`).
#;-) Warning: Removed 5 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/W9Oy1qP.png)<!-- -->

``` r
gof_res$Fetal_gof_res %>%
  ggplot(aes(Stage, r, fill = Stage)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples, correlation with Stage") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 5 rows containing non-finite values (`stat_ydensity()`).
#;-) Warning: Removed 5 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/EBcvmoR.png)<!-- -->

``` r
gof_res$Fetal_gof_res %>%
  ggplot(aes(Stage, r, fill = Regions)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  ggtitle("Goodness of fit for prenatal samples, correlation with regions and stage") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 5 rows containing non-finite values (`stat_boxplot()`).
#;-) Warning: Removed 5 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/TwLiv44.png)<!-- -->

``` r
head(gof_res$Fetal_gof_res)
#;-)                  SampleID       rho         r      mae     rmse        nmae
#;-) 1     13058_8_pcw_Ocx_s2a 0.4561448 0.4509607 16.40883 67.04400 0.005556920
#;-) 2 13058_8_pcw_M1C-S1C_s2a 0.5392566 0.5613217 15.04895 64.74304 0.005095913
#;-) 3     13058_8_pcw_AMY_s2a 0.4827804 0.4792736 16.04927 65.57019 0.005435154
#;-) 4     13058_8_pcw_MGE_s2a 0.5664574 0.5851683 15.16437 68.41082 0.005134971
#;-) 5     13058_8_pcw_STC_s2a 0.5469196 0.5647466 15.29077 68.93227 0.005177775
#;-) 6     13058_8_pcw_URL_s2a 0.5375673 0.5551688 15.46557 69.29532 0.005236966
#;-)            sigs X.1 column_num DonorID   DonorName   Age Sex StructureID
#;-) 1 Fetal_PFC_sig   1          1   13058 H376_IIA_51 8_pcw   M       10268
#;-) 2 Fetal_PFC_sig   2          2   13058 H376_IIA_51 8_pcw   M       10291
#;-) 3 Fetal_PFC_sig   3          3   13058 H376_IIA_51 8_pcw   M       10361
#;-) 4 Fetal_PFC_sig   4          4   13058 H376_IIA_51 8_pcw   M       10550
#;-) 5 Fetal_PFC_sig   5          5   13058 H376_IIA_51 8_pcw   M       10243
#;-) 6 Fetal_PFC_sig   6          6   13058 H376_IIA_51 8_pcw   M       10665
#;-)   StructureAcronym                                              Structure Stage
#;-) 1              Ocx                                    occipital neocortex   s2a
#;-) 2          M1C-S1C                 primary motor-sensory cortex (samples)   s2a
#;-) 3              AMY                                     amygdaloid complex   s2a
#;-) 4              MGE                             medial ganglionic eminence   s2a
#;-) 5              STC posterior (caudal) superior temporal cortex (area 22c)   s2a
#;-) 6              URL                            upper (rostral) rhombic lip   s2a
#;-)      Regions AgeInterval Diagnosis age_for_mRIN                   sample.name
#;-) 1     Cortex      8-9pcw   Control         8pcw     H376_IIA_51//8pcw//M//Ocx
#;-) 2     Cortex      8-9pcw   Control         8pcw H376_IIA_51//8pcw//M//M1C-S1C
#;-) 3  Subcortex      8-9pcw   Control         8pcw     H376_IIA_51//8pcw//M//AMY
#;-) 4  Subcortex      8-9pcw   Control         8pcw     H376_IIA_51//8pcw//M//MGE
#;-) 5     Cortex      8-9pcw   Control         8pcw     H376_IIA_51//8pcw//M//STC
#;-) 6 Cerebellum      8-9pcw   Control         8pcw     H376_IIA_51//8pcw//M//URL
#;-)         mRIN   z.score   P.value AgeNumeric   Period
#;-) 1 0.01961997 0.4167733 0.6626759 -0.6153846 Prenatal
#;-) 2         NA        NA        NA -0.6153846 Prenatal
#;-) 3 0.04569583 1.2579299 0.8967126 -0.6153846 Prenatal
#;-) 4 0.01104168 0.1400543 0.5566043 -0.6153846 Prenatal
#;-) 5 0.02622386 0.6298019 0.7367440 -0.6153846 Prenatal
#;-) 6 0.03846882 1.0248007 0.8483370 -0.6153846 Prenatal
```

``` r
gof_res$Adult_gof_res %>%
  rbind(., gof_res$Adult_gof_MB) %>%
  ggplot(aes(sigs, r, fill = sigs)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for postnatal samples") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
```

![](https://i.imgur.com/ClE2O7v.png)<!-- -->

``` r

gof_res$Adult_gof_MB %>%
  ggplot(aes(Regions, r, fill = Regions)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples, correlation with regions and stage") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
```

![](https://i.imgur.com/kXo0Rl0.png)<!-- -->

``` r

gof_res$Adult_gof_MB$Stage <- factor(gof_res$Adult_gof_MB$Stage, levels = c(
  "s6", "s7",
  "s8", "s9", "s10",
  "s11", "s12", "s13"
))



gof_res$Adult_gof_MB %>%
  ggplot(aes(Stage, r, fill = Regions)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  ggtitle("Goodness of fit for prenatal samples, correlation with regions and stage") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
```

![](https://i.imgur.com/3tB6kRN.png)<!-- -->

``` r
stage.md <- bspan.md %>% dplyr::filter(Period == "Postnatal")

max_est_adult <- apply(dtg_res$Adult[, -6], 1, max)

max_est_adult ## dtanlge not doing well
#;-)      13058_8_pcw_URL_s2a      13058_8_pcw_CGE_s2a      13058_8_pcw_DTH_s2a 
#;-)                0.9998924                0.9983001                0.9981654 
#;-)      13058_8_pcw_MFC_s2a      13058_8_pcw_DFC_s2a      13058_8_pcw_OFC_s2a 
#;-)                0.9999973                1.0000000                0.9999807 
#;-)      13058_8_pcw_LGE_s2a      13058_8_pcw_ITC_s2a      13058_8_pcw_HIP_s2a 
#;-)                0.9208008                0.9997899                0.9999810 
#;-)      13058_8_pcw_VFC_s2a      13058_8_pcw_PCx_s2a      12833_9_pcw_DFC_s2a 
#;-)                0.8965128                0.9990879                0.9999999 
#;-)      12833_9_pcw_MFC_s2a      12833_9_pcw_AMY_s2a      12833_9_pcw_DTH_s2a 
#;-)                0.9999993                0.8665887                0.9999996 
#;-)      12833_9_pcw_URL_s2a      12833_9_pcw_LGE_s2a  12833_9_pcw_M1C-S1C_s2a 
#;-)                0.9999660                0.9922629                0.9999999 
#;-)      12833_9_pcw_MGE_s2a      12833_9_pcw_TCx_s2a      12833_9_pcw_HIP_s2a 
#;-)                0.9999983                0.9988405                0.9959757 
#;-)      12833_9_pcw_Ocx_s2a      12833_9_pcw_CGE_s2a      12833_9_pcw_OFC_s2a 
#;-)                0.9999999                0.9999986                0.9999999 
#;-)      12833_9_pcw_PCx_s2a     12835_12_pcw_HIP_s2b     12835_12_pcw_DTH_s2b 
#;-)                0.6410536                0.9998106                1.0000000 
#;-)     12835_12_pcw_AMY_s2b     12835_12_pcw_A1C_s2b     12835_12_pcw_V1C_s2b 
#;-)                0.8266168                1.0000000                1.0000000 
#;-)     12835_12_pcw_MFC_s2b     12835_12_pcw_STR_s2b     12835_12_pcw_M1C_s2b 
#;-)                0.9999765                1.0000000                1.0000000 
#;-)     12835_12_pcw_DFC_s2b     12835_12_pcw_ITC_s2b     12835_12_pcw_VFC_s2b 
#;-)                1.0000000                0.9951060                1.0000000 
#;-)     12835_12_pcw_STC_s2b     12835_12_pcw_IPC_s2b     12835_12_pcw_OFC_s2b 
#;-)                1.0000000                1.0000000                0.9999997 
#;-)     12835_12_pcw_S1C_s2b     12960_12_pcw_STC_s2b     12960_12_pcw_V1C_s2b 
#;-)                0.9999999                1.0000000                1.0000000 
#;-)     12960_12_pcw_A1C_s2b      12960_12_pcw_CB_s2b     12960_12_pcw_ITC_s2b 
#;-)                1.0000000                0.9998696                1.0000000 
#;-)     12960_12_pcw_IPC_s2b     12960_12_pcw_S1C_s2b     12960_12_pcw_DTH_s2b 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)     12960_12_pcw_HIP_s2b     12960_12_pcw_DFC_s2b     12960_12_pcw_OFC_s2b 
#;-)                0.9400885                1.0000000                0.9997545 
#;-)     12960_12_pcw_M1C_s2b     12960_12_pcw_VFC_s2b     12960_12_pcw_AMY_s2b 
#;-)                1.0000000                0.9999986                0.9999983 
#;-)     12960_12_pcw_STR_s2b     13060_12_pcw_AMY_s2b     13060_12_pcw_M1C_s2b 
#;-)                1.0000000                1.0000000                0.9999991 
#;-)     13060_12_pcw_ITC_s2b     13060_12_pcw_VFC_s2b     13060_12_pcw_STR_s2b 
#;-)                0.9819100                0.9979371                1.0000000 
#;-)     13060_12_pcw_CBC_s2b     13060_12_pcw_DTH_s2b     13060_12_pcw_DFC_s2b 
#;-)                1.0000000                1.0000000                0.9999956 
#;-)     13060_12_pcw_OFC_s2b     13060_12_pcw_S1C_s2b     13060_12_pcw_V1C_s2b 
#;-)                0.9999990                0.9999686                0.9999999 
#;-)     13060_12_pcw_MFC_s2b     13060_12_pcw_HIP_s2b     13060_12_pcw_IPC_s2b 
#;-)                0.9999991                0.6318030                0.9796768 
#;-)     13060_12_pcw_A1C_s2b     12820_13_pcw_STR_s3a     12820_13_pcw_M1C_s3a 
#;-)                0.9995040                0.9999996                1.0000000 
#;-)     12820_13_pcw_AMY_s3a     12820_13_pcw_S1C_s3a     12820_13_pcw_V1C_s3a 
#;-)                0.9991861                0.9999999                1.0000000 
#;-)     12820_13_pcw_A1C_s3a     12820_13_pcw_HIP_s3a     12820_13_pcw_VFC_s3a 
#;-)                0.9992115                0.9043901                1.0000000 
#;-)     12820_13_pcw_MFC_s3a     12820_13_pcw_STC_s3a     12820_13_pcw_OFC_s3a 
#;-)                0.9999835                1.0000000                0.9999828 
#;-)     12820_13_pcw_IPC_s3a     12820_13_pcw_DFC_s3a     12820_13_pcw_ITC_s3a 
#;-)                1.0000000                1.0000000                0.9999988 
#;-)     12834_13_pcw_STR_s3a     12834_13_pcw_S1C_s3a     12834_13_pcw_V1C_s3a 
#;-)                1.0000000                0.9999998                0.9999998 
#;-)     12834_13_pcw_AMY_s3a     12834_13_pcw_IPC_s3a      12834_13_pcw_CB_s3a 
#;-)                0.9995919                0.9998592                0.9999537 
#;-)     12834_13_pcw_DFC_s3a     12834_13_pcw_OFC_s3a     12834_13_pcw_STC_s3a 
#;-)                0.9999550                0.9999941                0.9834309 
#;-)     12834_13_pcw_A1C_s3a     12834_13_pcw_HIP_s3a     12834_13_pcw_M1C_s3a 
#;-)                0.9510396                0.9264854                0.9999993 
#;-)     12834_13_pcw_ITC_s3a      12834_13_pcw_MD_s3a     12834_13_pcw_MFC_s3a 
#;-)                0.9004478                1.0000000                0.9999932 
#;-)     12834_13_pcw_VFC_s3a     12888_13_pcw_S1C_s3a     12888_13_pcw_DFC_s3a 
#;-)                0.9999996                0.9999858                0.9990413 
#;-)      12888_13_pcw_CB_s3a     12888_13_pcw_AMY_s3a     12888_13_pcw_ITC_s3a 
#;-)                0.9999667                0.9922512                0.9999362 
#;-)     12888_13_pcw_HIP_s3a     12888_13_pcw_IPC_s3a     12888_13_pcw_A1C_s3a 
#;-)                0.8298414                0.9999999                0.9299502 
#;-)     12888_13_pcw_VFC_s3a     12888_13_pcw_STR_s3a     12888_13_pcw_OFC_s3a 
#;-)                0.9975760                1.0000000                0.8161353 
#;-)     12888_13_pcw_MFC_s3a     12888_13_pcw_V1C_s3a     12888_13_pcw_M1C_s3a 
#;-)                0.7902382                0.9997887                0.9999952 
#;-)     12287_16_pcw_VFC_s3b     12287_16_pcw_OFC_s3b     12287_16_pcw_STC_s3b 
#;-)                0.9968546                0.9773448                0.9997985 
#;-)     12287_16_pcw_DFC_s3b     12287_16_pcw_IPC_s3b     12287_16_pcw_A1C_s3b 
#;-)                0.9991808                0.8871932                0.5831539 
#;-)     12287_16_pcw_MFC_s3b     12287_16_pcw_STR_s3b      12287_16_pcw_MD_s3b 
#;-)                0.9991406                1.0000000                1.0000000 
#;-)     12287_16_pcw_V1C_s3b     12837_16_pcw_VFC_s3b     12837_16_pcw_STC_s3b 
#;-)                0.8980962                1.0000000                0.9999995 
#;-)     12837_16_pcw_IPC_s3b     12837_16_pcw_A1C_s3b      12837_16_pcw_MD_s3b 
#;-)                0.9999997                0.9999998                0.9999982 
#;-)     12837_16_pcw_AMY_s3b     12837_16_pcw_V1C_s3b     12837_16_pcw_ITC_s3b 
#;-)                0.7303595                0.9999982                0.9993695 
#;-)     12837_16_pcw_CBC_s3b     12837_16_pcw_HIP_s3b     12837_16_pcw_STR_s3b 
#;-)                0.9771349                0.9901485                0.9999994 
#;-)     12837_16_pcw_OFC_s3b     12837_16_pcw_DFC_s3b     12837_16_pcw_M1C_s3b 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)     12837_16_pcw_MFC_s3b     12837_16_pcw_S1C_s3b     12879_16_pcw_HIP_s3b 
#;-)                1.0000000                0.9999999                0.9876680 
#;-)     12879_16_pcw_ITC_s3b     12879_16_pcw_DFC_s3b     12879_16_pcw_IPC_s3b 
#;-)                0.9986692                0.9355196                0.9999945 
#;-)     12879_16_pcw_STC_s3b     12879_16_pcw_V1C_s3b      12879_16_pcw_MD_s3b 
#;-)                0.9820298                0.9319316                1.0000000 
#;-)     12879_16_pcw_STR_s3b 12879_16_pcw_M1C-S1C_s3b     12879_16_pcw_A1C_s3b 
#;-)                1.0000000                0.9999234                0.5863975 
#;-)     12879_16_pcw_VFC_s3b     12879_16_pcw_AMY_s3b     12879_16_pcw_MFC_s3b 
#;-)                0.9999930                0.9981499                0.9999988 
#;-)     12880_17_pcw_STR_s3b     12880_17_pcw_STC_s3b      12880_17_pcw_MD_s3b 
#;-)                1.0000000                0.5332172                1.0000000 
#;-)     12880_17_pcw_A1C_s3b     12880_17_pcw_HIP_s3b     12880_17_pcw_DFC_s3b 
#;-)                0.8569374                0.9981833                0.8307139 
#;-)     12880_17_pcw_VFC_s3b     12880_17_pcw_IPC_s3b 12880_17_pcw_M1C-S1C_s3b 
#;-)                0.9528005                0.9826520                0.9654870 
#;-)     12880_17_pcw_V1C_s3b     12880_17_pcw_CBC_s3b     12880_17_pcw_OFC_s3b 
#;-)                0.9547846                0.8367436                0.7602700 
#;-)     12880_17_pcw_AMY_s3b     12880_17_pcw_MFC_s3b      12885_19_pcw_STC_s4 
#;-)                0.4610854                0.9400968                0.9999952 
#;-)      12885_19_pcw_IPC_s4  12885_19_pcw_M1C-S1C_s4      12885_19_pcw_MFC_s4 
#;-)                0.9999983                0.9999991                0.9999280 
#;-)      12885_19_pcw_V1C_s4      12885_19_pcw_HIP_s4      12885_19_pcw_A1C_s4 
#;-)                0.9999954                0.9997490                0.9999883 
#;-)       12885_19_pcw_MD_s4      12885_19_pcw_VFC_s4      12885_19_pcw_STR_s4 
#;-)                1.0000000                0.9894958                0.9999999 
#;-)      12885_19_pcw_DFC_s4      12365_21_pcw_CBC_s4      12365_21_pcw_ITC_s4 
#;-)                0.9999345                0.8491391                0.9973377 
#;-)      12886_21_pcw_MFC_s4      12886_21_pcw_M1C_s4      12886_21_pcw_STR_s4 
#;-)                0.6492185                0.7012891                0.9999997 
#;-)      12886_21_pcw_HIP_s4      12886_21_pcw_STC_s4      12886_21_pcw_ITC_s4 
#;-)                0.9213276                0.9886051                0.6253089 
#;-)      12886_21_pcw_AMY_s4      12886_21_pcw_VFC_s4      12886_21_pcw_OFC_s4 
#;-)                0.9672540                0.9034668                0.9818229 
#;-)      12886_21_pcw_V1C_s4      12886_21_pcw_S1C_s4      12886_21_pcw_CBC_s4 
#;-)                0.9585501                0.6175358                0.9666068 
#;-)      12886_21_pcw_IPC_s4      12886_21_pcw_DFC_s4      12288_24_pcw_HIP_s4 
#;-)                0.9526652                0.9625317                0.9996480 
#;-)      12288_24_pcw_AMY_s4      12288_24_pcw_DFC_s4      12288_24_pcw_S1C_s4 
#;-)                0.9880562                0.9668620                0.6404061 
#;-)       12288_24_pcw_MD_s4      12288_24_pcw_STR_s4      12288_24_pcw_A1C_s4 
#;-)                1.0000000                0.9999994                0.9586770 
#;-)      12288_24_pcw_IPC_s4      12288_24_pcw_OFC_s4      12288_24_pcw_STC_s4 
#;-)                0.9972587                0.9991387                0.9958906 
#;-)      12288_24_pcw_ITC_s4      12288_24_pcw_VFC_s4      12288_24_pcw_M1C_s4 
#;-)                0.9561578                0.9861831                0.9815373 
#;-)      12288_24_pcw_MFC_s4      12288_24_pcw_CBC_s4      12288_24_pcw_V1C_s4 
#;-)                0.9918946                0.9920447                0.6567063 
#;-)      12948_25_pcw_A1C_s5      12949_26_pcw_STC_s5      12949_26_pcw_V1C_s5 
#;-)                1.0000000                0.6317932                1.0000000 
#;-)      12949_26_pcw_DFC_s5      12295_35_pcw_CBC_s5      12295_35_pcw_VFC_s5 
#;-)                1.0000000                0.6894398                0.9999990 
#;-)  263195015_37_pcw_A1C_s5  263195015_37_pcw_OFC_s5  263195015_37_pcw_AMY_s5 
#;-)                0.9998357                0.9888975                0.9999816 
#;-)  263195015_37_pcw_V1C_s5  263195015_37_pcw_MFC_s5   263195015_37_pcw_MD_s5 
#;-)                0.9997018                0.9996289                0.9999999 
#;-)  263195015_37_pcw_CBC_s5  263195015_37_pcw_STR_s5  263195015_37_pcw_IPC_s5 
#;-)                1.0000000                1.0000000                0.9987284 
#;-)  263195015_37_pcw_ITC_s5  263195015_37_pcw_S1C_s5  263195015_37_pcw_DFC_s5 
#;-)                0.9913987                0.9999538                0.9741979 
#;-)  263195015_37_pcw_STC_s5  263195015_37_pcw_M1C_s5  263195015_37_pcw_VFC_s5 
#;-)                0.9899308                0.9999999                0.9999015 
#;-)  263195015_37_pcw_HIP_s5       12296_4_mos_STC_s6       12296_4_mos_STR_s6 
#;-)                1.0000000                0.9747563                0.9490654 
#;-)       12296_4_mos_V1C_s6       12296_4_mos_CBC_s6       12296_4_mos_HIP_s6 
#;-)                0.7223764                0.9716594                0.7047771 
#;-)        12296_4_mos_MD_s6       12296_4_mos_ITC_s6       12296_4_mos_AMY_s6 
#;-)                0.9999988                0.9986360                0.5779782 
#;-)       12296_4_mos_MFC_s6       12889_4_mos_M1C_s6       12889_4_mos_DFC_s6 
#;-)                0.9998039                1.0000000                0.5764035 
#;-)       12889_4_mos_OFC_s6       12889_4_mos_A1C_s6       12889_4_mos_STC_s6 
#;-)                0.9999996                0.9976647                0.9693527 
#;-)       12889_4_mos_VFC_s6       12889_4_mos_ITC_s6       12889_4_mos_AMY_s6 
#;-)                0.9966422                0.9935795                0.9999991 
#;-)       12890_4_mos_M1C_s6        12890_4_mos_MD_s6       12890_4_mos_OFC_s6 
#;-)                0.8723398                0.9996106                0.9999921 
#;-)       12890_4_mos_STR_s6       12890_4_mos_STC_s6       12890_4_mos_A1C_s6 
#;-)                0.9999998                0.9720623                0.7430429 
#;-)       12890_4_mos_HIP_s6       12890_4_mos_AMY_s6       12890_4_mos_MFC_s6 
#;-)                0.9999999                0.9981220                0.9797703 
#;-)       12890_4_mos_V1C_s6       12890_4_mos_ITC_s6       12890_4_mos_S1C_s6 
#;-)                0.5045675                0.9980102                0.9940263 
#;-)       12890_4_mos_IPC_s6       12890_4_mos_CBC_s6       12890_4_mos_DFC_s6 
#;-)                0.9592347                0.9999949                0.9772912 
#;-)       12890_4_mos_VFC_s6      12977_10_mos_S1C_s7      12977_10_mos_IPC_s7 
#;-)                0.9530646                0.9997900                0.9999347 
#;-)      12977_10_mos_STC_s7      12977_10_mos_DFC_s7      12977_10_mos_OFC_s7 
#;-)                0.9903027                0.9999896                0.9925425 
#;-)      12977_10_mos_MFC_s7      12977_10_mos_CBC_s7      12977_10_mos_ITC_s7 
#;-)                0.8659817                0.9973495                0.9859711 
#;-)       12977_10_mos_MD_s7      12977_10_mos_V1C_s7       12830_1_yrs_AMY_s7 
#;-)                0.9982553                0.9764874                0.9942854 
#;-)       12830_1_yrs_M1C_s7       12830_1_yrs_DFC_s7       12830_1_yrs_STC_s7 
#;-)                0.9999810                0.9957497                0.9998584 
#;-)       12830_1_yrs_V1C_s7       12830_1_yrs_MFC_s7       12830_1_yrs_OFC_s7 
#;-)                0.9999600                0.9966607                0.7301544 
#;-)       12830_1_yrs_ITC_s7       12830_1_yrs_S1C_s7       12830_1_yrs_VFC_s7 
#;-)                0.9975058                0.9999643                0.9999458 
#;-)       12830_1_yrs_CBC_s7       12830_1_yrs_A1C_s7        12830_1_yrs_MD_s7 
#;-)                0.9999936                0.9999251                0.9988382 
#;-)       12830_1_yrs_HIP_s7       12830_1_yrs_IPC_s7       12830_1_yrs_STR_s7 
#;-)                0.9796795                0.9995016                0.9991174 
#;-)       12979_2_yrs_VFC_s8       12979_2_yrs_IPC_s8       12979_2_yrs_S1C_s8 
#;-)                0.9999945                1.0000000                1.0000000 
#;-)       12979_2_yrs_ITC_s8        12979_2_yrs_MD_s8       12979_2_yrs_DFC_s8 
#;-)                1.0000000                1.0000000                0.9999994 
#;-)       12979_2_yrs_MFC_s8       12979_2_yrs_CBC_s8       12979_2_yrs_HIP_s8 
#;-)                1.0000000                1.0000000                0.9997881 
#;-)       12979_2_yrs_OFC_s8       12979_2_yrs_V1C_s8       12979_2_yrs_STC_s8 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)       12836_3_yrs_STR_s8       12836_3_yrs_STC_s8       12836_3_yrs_IPC_s8 
#;-)                0.9999948                0.9950686                0.6706557 
#;-)       12836_3_yrs_ITC_s8       12836_3_yrs_CBC_s8       12836_3_yrs_AMY_s8 
#;-)                0.9998517                1.0000000                0.9981214 
#;-)       12836_3_yrs_VFC_s8       12836_3_yrs_M1C_s8       12836_3_yrs_A1C_s8 
#;-)                0.9601390                0.9919670                0.9067903 
#;-)        12836_3_yrs_MD_s8       12836_3_yrs_V1C_s8       12980_3_yrs_OFC_s8 
#;-)                0.9998994                0.9979865                1.0000000 
#;-)       12980_3_yrs_MFC_s8       12980_3_yrs_M1C_s8       12980_3_yrs_S1C_s8 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)       12980_3_yrs_STC_s8       12980_3_yrs_IPC_s8       12980_3_yrs_A1C_s8 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)       12980_3_yrs_CBC_s8       12980_3_yrs_V1C_s8       12980_3_yrs_AMY_s8 
#;-)                1.0000000                1.0000000                1.0000000 
#;-)       12980_3_yrs_VFC_s8       12980_3_yrs_DFC_s8       12980_3_yrs_HIP_s8 
#;-)                1.0000000                1.0000000                0.9999999 
#;-)       12980_3_yrs_ITC_s8       12298_4_yrs_AMY_s8       12298_4_yrs_STR_s8 
#;-)                1.0000000                0.9999999                0.9962710 
#;-)       12298_4_yrs_VFC_s8       12298_4_yrs_DFC_s8       12298_4_yrs_STC_s8 
#;-)                0.9999994                1.0000000                0.9999999 
#;-)       12298_4_yrs_CBC_s8        12298_4_yrs_MD_s8       12841_8_yrs_ITC_s9 
#;-)                0.9937467                0.9999882                0.9995397 
#;-)       12841_8_yrs_S1C_s9       12841_8_yrs_IPC_s9       12841_8_yrs_V1C_s9 
#;-)                0.9999989                0.9999587                0.9988656 
#;-)       12841_8_yrs_STC_s9        12841_8_yrs_MD_s9       12841_8_yrs_VFC_s9 
#;-)                0.9999293                0.9999982                0.9999943 
#;-)       12841_8_yrs_HIP_s9       12841_8_yrs_AMY_s9       12841_8_yrs_M1C_s9 
#;-)                0.9984386                0.8258419                0.9999975 
#;-)       12841_8_yrs_CBC_s9       12841_8_yrs_DFC_s9       12841_8_yrs_STR_s9 
#;-)                0.9783105                0.9962564                0.9999671 
#;-)       12841_8_yrs_A1C_s9       12841_8_yrs_MFC_s9       12841_8_yrs_OFC_s9 
#;-)                0.9998660                0.9999178                0.9717242 
#;-)       12981_8_yrs_MFC_s9       12981_8_yrs_CBC_s9       12981_8_yrs_A1C_s9 
#;-)                0.9991644                1.0000000                0.9999996 
#;-)       12981_8_yrs_VFC_s9       12981_8_yrs_HIP_s9       12981_8_yrs_AMY_s9 
#;-)                0.9902094                0.9985784                0.9999954 
#;-)       12981_8_yrs_V1C_s9       12981_8_yrs_IPC_s9       12981_8_yrs_DFC_s9 
#;-)                0.9999537                0.9999995                0.9967980 
#;-)       12981_8_yrs_STC_s9       12981_8_yrs_ITC_s9      12289_11_yrs_VFC_s9 
#;-)                0.9999924                1.0000000                0.9657305 
#;-)      12289_11_yrs_IPC_s9      12289_11_yrs_ITC_s9      12289_11_yrs_STC_s9 
#;-)                0.9833616                0.9827440                0.9876206 
#;-)      12289_11_yrs_S1C_s9      12289_11_yrs_DFC_s9      12289_11_yrs_AMY_s9 
#;-)                0.5531076                0.7436341                0.9935592 
#;-)      12289_11_yrs_V1C_s9      12289_11_yrs_CBC_s9      12289_11_yrs_A1C_s9 
#;-)                0.9941382                0.9999999                0.9998976 
#;-)      12289_11_yrs_HIP_s9      12289_11_yrs_M1C_s9      12289_11_yrs_OFC_s9 
#;-)                0.9999999                0.6645086                0.9294756 
#;-)      12289_11_yrs_MFC_s9     12831_13_yrs_MFC_s10     12831_13_yrs_IPC_s10 
#;-)                0.7611074                0.9998334                0.9999792 
#;-)     12831_13_yrs_VFC_s10     12831_13_yrs_A1C_s10     12831_13_yrs_ITC_s10 
#;-)                0.9997659                0.9999888                0.9995243 
#;-)     12831_13_yrs_STC_s10     12831_13_yrs_CBC_s10     12831_13_yrs_V1C_s10 
#;-)                0.9782420                1.0000000                0.9998945 
#;-)     12831_13_yrs_S1C_s10     12831_13_yrs_DFC_s10     12831_13_yrs_OFC_s10 
#;-)                0.9990813                0.9999940                0.9999683 
#;-)      12831_13_yrs_MD_s10     12831_13_yrs_AMY_s10     12831_13_yrs_M1C_s10 
#;-)                0.9999998                0.9996745                0.9998261 
#;-)     12831_13_yrs_HIP_s10     12831_13_yrs_STR_s10     12299_15_yrs_IPC_s10 
#;-)                0.9999996                0.9999984                0.9996853 
#;-)     12299_15_yrs_CBC_s10     12299_15_yrs_STC_s10     12299_15_yrs_ITC_s10 
#;-)                0.9999792                0.9677742                0.9844680 
#;-)     12299_15_yrs_AMY_s10     12984_18_yrs_HIP_s10     12984_18_yrs_M1C_s10 
#;-)                0.9998990                0.9999900                0.9999925 
#;-)     12984_18_yrs_V1C_s10     12984_18_yrs_A1C_s10     12984_18_yrs_ITC_s10 
#;-)                0.9999998                0.9999986                0.9999189 
#;-)     12984_18_yrs_VFC_s10     12984_18_yrs_MFC_s10     12984_18_yrs_STC_s10 
#;-)                0.9998616                0.9977959                0.9999828 
#;-)     12984_18_yrs_DFC_s10     12984_18_yrs_S1C_s10     12984_18_yrs_IPC_s10 
#;-)                0.9999096                0.9999853                0.9993271 
#;-)     12984_18_yrs_OFC_s10     12984_18_yrs_CBC_s10     12832_19_yrs_STC_s10 
#;-)                0.9962075                0.9999962                0.9583503 
#;-)     12832_19_yrs_CBC_s10     12832_19_yrs_DFC_s10     12832_19_yrs_AMY_s10 
#;-)                0.9991536                0.7496035                0.9689095 
#;-)     12832_19_yrs_OFC_s10     12832_19_yrs_ITC_s10     12832_19_yrs_STR_s10 
#;-)                0.5469061                0.9512882                0.9999487 
#;-)     12832_19_yrs_S1C_s10     12832_19_yrs_M1C_s10     12832_19_yrs_HIP_s10 
#;-)                0.9568104                0.7275802                0.8714747 
#;-)     12832_19_yrs_A1C_s10     12832_19_yrs_V1C_s10     12832_19_yrs_IPC_s10 
#;-)                0.9188647                0.6459457                0.9984479 
#;-)     12832_19_yrs_VFC_s10     12832_19_yrs_MFC_s10      12832_19_yrs_MD_s10 
#;-)                0.9885291                0.9860401                0.9999895 
#;-)     13057_21_yrs_MFC_s11     13057_21_yrs_M1C_s11     13057_21_yrs_A1C_s11 
#;-)                0.9999936                0.9999828                0.9999003 
#;-)     13057_21_yrs_DFC_s11     13057_21_yrs_AMY_s11      13057_21_yrs_MD_s11 
#;-)                0.9999995                0.9974679                0.9596877 
#;-)     13057_21_yrs_V1C_s11     13057_21_yrs_STC_s11     13057_21_yrs_VFC_s11 
#;-)                0.9999262                0.9999698                0.9999947 
#;-)     13057_21_yrs_CBC_s11     13057_21_yrs_S1C_s11     13057_21_yrs_ITC_s11 
#;-)                1.0000000                0.9999977                0.9975536 
#;-)     13057_21_yrs_IPC_s11     13057_21_yrs_HIP_s11     13057_21_yrs_STR_s11 
#;-)                0.9999990                0.9851133                0.9999902 
#;-)     13057_21_yrs_OFC_s11     12300_23_yrs_ITC_s11     12300_23_yrs_AMY_s11 
#;-)                0.9999744                0.9999195                0.9999999 
#;-)     12300_23_yrs_OFC_s11      12300_23_yrs_MD_s11     12300_23_yrs_STC_s11 
#;-)                0.9999925                1.0000000                0.9999983 
#;-)     12300_23_yrs_VFC_s11     12300_23_yrs_MFC_s11     12300_23_yrs_M1C_s11 
#;-)                0.9996442                0.9999307                0.9998713 
#;-)     12300_23_yrs_HIP_s11     12300_23_yrs_CBC_s11     12300_23_yrs_S1C_s11 
#;-)                0.9999988                0.9999999                0.9998206 
#;-)     12300_23_yrs_A1C_s11     12300_23_yrs_IPC_s11     12300_23_yrs_STR_s11 
#;-)                0.9999919                0.9999406                0.9999799 
#;-)     12290_30_yrs_ITC_s12     12290_30_yrs_STR_s12     12290_30_yrs_MFC_s12 
#;-)                0.9628903                0.9292228                0.6542077 
#;-)     12290_30_yrs_OFC_s12     12290_30_yrs_IPC_s12     12290_30_yrs_STC_s12 
#;-)                0.9566861                0.7468012                0.5191693 
#;-)     12290_30_yrs_S1C_s12      12290_30_yrs_MD_s12     12290_30_yrs_HIP_s12 
#;-)                0.7187381                0.9999987                0.9731308 
#;-)     12290_30_yrs_A1C_s12     12290_30_yrs_V1C_s12     12290_30_yrs_DFC_s12 
#;-)                0.7535591                0.9688586                0.9876675 
#;-)     12290_30_yrs_AMY_s12     12290_30_yrs_M1C_s12     12290_30_yrs_VFC_s12 
#;-)                0.9732397                0.7082106                0.6305485 
#;-)     12290_30_yrs_CBC_s12     12302_36_yrs_M1C_s12     12302_36_yrs_STC_s12 
#;-)                0.9734140                0.9928672                0.9150881 
#;-)     12302_36_yrs_MFC_s12     12302_36_yrs_ITC_s12     12302_36_yrs_IPC_s12 
#;-)                0.9999831                0.9763748                0.9814993 
#;-)     12302_36_yrs_AMY_s12     12302_36_yrs_OFC_s12     12302_36_yrs_VFC_s12 
#;-)                0.8734098                0.9177011                0.7976362 
#;-)     12302_36_yrs_DFC_s12     12302_36_yrs_A1C_s12     12302_36_yrs_V1C_s12 
#;-)                0.7882808                0.5146300                0.7412773 
#;-)     12302_36_yrs_CBC_s12     12302_36_yrs_S1C_s12     12302_36_yrs_HIP_s12 
#;-)                0.9999840                0.9924636                0.8883397 
#;-)     12302_36_yrs_STR_s12      12302_36_yrs_MD_s12     12303_37_yrs_S1C_s12 
#;-)                0.9582242                0.9999821                0.9699596 
#;-)     12303_37_yrs_V1C_s12     12303_37_yrs_CBC_s12     12303_37_yrs_HIP_s12 
#;-)                0.9789559                0.6776341                0.9733735 
#;-)     12303_37_yrs_IPC_s12      12303_37_yrs_MD_s12     12303_37_yrs_DFC_s12 
#;-)                0.9489543                0.9999973                0.9953689 
#;-)     12303_37_yrs_VFC_s12     12303_37_yrs_ITC_s12     12303_37_yrs_A1C_s12 
#;-)                0.9890883                0.9975659                0.9994768 
#;-)     12303_37_yrs_M1C_s12     12303_37_yrs_OFC_s12     12303_37_yrs_STR_s12 
#;-)                0.9991508                0.9975483                0.9773402 
#;-)     12303_37_yrs_STC_s12     12303_37_yrs_AMY_s12     12303_37_yrs_MFC_s12 
#;-)                0.6626766                0.9924880                0.9999992 
#;-)     12304_40_yrs_DFC_s13     12304_40_yrs_ITC_s13     12304_40_yrs_VFC_s13 
#;-)                0.9999999                0.9192647                0.6440340 
#;-)      12304_40_yrs_MD_s13     12304_40_yrs_AMY_s13     12304_40_yrs_A1C_s13 
#;-)                0.9868781                0.9999999                0.9999998 
#;-)     12304_40_yrs_CBC_s13     12304_40_yrs_V1C_s13     12304_40_yrs_OFC_s13 
#;-)                1.0000000                0.8047138                0.6313344 
#;-)     12304_40_yrs_STC_s13     12304_40_yrs_IPC_s13     12304_40_yrs_M1C_s13 
#;-)                0.5495787                0.6689818                1.0000000 
#;-)     12304_40_yrs_HIP_s13     12304_40_yrs_STR_s13     12304_40_yrs_S1C_s13 
#;-)                1.0000000                0.8110735                0.9999970 
#;-)              Astro.Adult              Micro.Adult              Oligo.Adult 
#;-)                0.8391702                0.9947247                0.9806886 
#;-)                OPC.Adult             Neuron.Adult 
#;-)                0.9980198                0.8381562

max_est_fetal %<>% melt() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% stage.md$SampleID) %>%
  left_join(gof_res$Fetal_gof_res)
#;-) Warning in melt(.): The melt generic in data.table has been passed a numeric
#;-) and will attempt to redirect to the relevant reshape2 method; please note that
#;-) reshape2 is deprecated, and this redirection is now deprecated as well. To
#;-) continue using melt methods from reshape2 while both libraries are attached,
#;-) e.g. melt.list, you can prepend the namespace like reshape2::melt(.). In the
#;-) next version, this warning will become an error.
#;-) Joining, by = "SampleID"

m <- match(gof_res$Fetal_gof_res$SampleID, names(max_est_fetal))

cor(max_est_fetal[m], gof_res$Fetal_gof_res$r, use = "pairwise.complete.obs")
#;-) Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'cor': undefined columns selected
```

## Rerunning this analysis, just with cortex samples

``` r
stages <- c("Fetal", "Adult")
gof_res_cortex <- list()

for (stage in stages) {
  if (stage == "Fetal") {
    stage.md <- bspan.md %>% dplyr::filter(Period == "Prenatal" & Regions == "Cortex")
  } else {
    stage.md <- bspan.md %>% dplyr::filter(Period == "Postnatal" & Regions == "Cortex")
  }

  stage.exp <- bspan.exp %>% dplyr::select(stage.md$SampleID)
  dtg.stage <- dtg_res[[paste0(stage)]] %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))
  gof <- write.gof(stage.exp, dtg.stage, sigs_dev) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_PFC_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)
  gof_res_cortex[[paste0(stage, "_gof_res")]] <- gof

  ## MultiBrain

  dtg.MB <- dtg_res$MB %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))

  gof.MB <- write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_MB_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)

  gof_res_cortex[[paste0(stage, "_gof_MB")]] <- gof.MB
}
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
```

``` r
gof_res_cortex$Fetal_gof_res %>%
  rbind(., gof_res_cortex$Fetal_gof_MB) %>%
  ggplot(aes(sigs, r, fill = sigs)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 6 rows containing non-finite values (`stat_ydensity()`).
#;-) Warning: Removed 6 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/hpmtoh2.png)<!-- -->

``` r
gof_res_cortex$Fetal_gof_res %>%
  ggplot(aes(StructureAcronym, r, fill = StructureAcronym)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
#;-) Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).
#;-) Warning: Groups with fewer than two data points have been dropped.
#;-) Warning: Removed 3 rows containing missing values (`geom_point()`).
```

![](https://i.imgur.com/RxcDI1o.png)<!-- -->

``` r
gof_res_cortex$Adult_gof_res %>%
  rbind(., gof_res_cortex$Adult_gof_MB) %>%
  ggplot(aes(sigs, r, fill = sigs)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Goodness of fit for prenatal samples") +
  ylab("Goodness of fit (r)") +
  xlab("Signatures")
```

![](https://i.imgur.com/HToN5O5.png)<!-- -->

- correlate gof with different brain regions?
- are the samples with a lower gof associated with different regions?
- adult brain samples? would it have better correlation?

# Deconvolution with dtangle

``` r
source("libraries.R")
source("functions.R")

load("../../Results/signatures/pfc_signatures.Rda")
load("../../DeconRNAShiny/sigsBrain.rda")
```

``` r
bspan.exp <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-exp.csv", row.names = 2, check.names = FALSE)[, -1]
bspan.md <- read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv", row.names = 1, header = TRUE)
```

``` r
stages.df <- list()
stages.df[["Fetal"]] <- c(paste(c(4:40), "pcw", sep = "_"))
stages.df[["Neonatal"]] <- c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Infancy"]] <- c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Childhood"]] <- paste(c(1:10), "yrs", sep = "_")
stages.df[["Adolescence"]] <- paste(c(11:19), "yrs", sep = "_")
stages.df[["Adult"]] <- paste(c(20:40), "yrs", sep = "_")


add_feature <- function(feature_column, features) {
  as.vector(sapply(feature_column, function(x) {
    names(features)[sapply(features, function(f) x %in% f)]
  }))
}

bspan.md <- bspan.md %>%
  mutate(Dev.stage = add_feature(.$Age, stages.df))

table(bspan.md$Dev.stage)
#;-) Error in base::table(...): all arguments must have the same length
```

``` r
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

dtg_res <- list()
for (stage in stages) {
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))

  dev_res <- run_dtg(bspan.exp, sigs_dev) %>%
    as.data.frame() %>%
    mutate(Sig = stage)

  dtg_res[[paste0(stage)]] <- dev_res
}
#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf




## 2) run dtangle for MB

dtg_res$MB <- run_dtg(bspan.exp, sigsBrain$MB) %>%
  as.data.frame() %>%
  mutate(Sig = "MB")
#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

#;-) Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf
```

``` r
stages <- c("Fetal", "Neonatal", "Childhood", "Adolescence", "Adult")
gof_res <- list()

for (stage in stages) {
  stage.md <- bspan.md %>% dplyr::filter(Dev.stage == stage)
  stage.exp <- bspan.exp %>% dplyr::select(stage.md$SampleID)
  dtg.stage <- dtg_res[[paste0(stage)]] %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))
  gof <- write.gof(stage.exp, dtg.stage, sigs_dev) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_PFC_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)
  gof_res[[paste0(stage, "_gof_res")]] <- gof

  ## MultiBrain

  dtg.MB <- dtg_res$MB %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))

  gof.MB <- write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_MB_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)

  gof_res[[paste0(stage, "_gof_MB")]] <- gof.MB
}
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': attempt to select less than one element in integerOneIndex
```

``` r
stages.df <- list()
stages.df[["Fetal"]] <- c(paste(c(4:40), "pcw", sep = "_"))
# stages.df[["Neonatal"]] =c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Infancy"]] <- c(paste(c(0:10), "mos", sep = "_"))
stages.df[["Childhood"]] <- paste(c(1:10), "yrs", sep = "_")
stages.df[["Adolescence"]] <- paste(c(11:19), "yrs", sep = "_")
stages.df[["Adult"]] <- paste(c(20:40), "yrs", sep = "_")


add_feature <- function(feature_column, features) {
  as.vector(sapply(feature_column, function(x) {
    names(features)[sapply(features, function(f) x %in% f)]
  }))
}

bspan.md <- bspan.md %>%
  mutate(Dev.stage = add_feature(.$Age, stages.df))

table(bspan.md$Dev.stage)
#;-) 
#;-) Adolescence       Adult   Childhood       Fetal     Infancy 
#;-)          64          93          87         237          43
```

``` r
stages <- c("Fetal", "Infancy", "Childhood", "Adolescence", "Adult")
for (stage in stages) {
  stage.md <- bspan.md %>% dplyr::filter(Dev.stage == stage)
  stage.exp <- bspan.exp %>% dplyr::select(stage.md$SampleID)
  dtg.stage <- dtg_res[[paste0(stage)]] %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))
  sigs_dev <- pfc_signatures$rpkm_all_neuro %>%
    dplyr::select(contains(stage))
  gof <- write.gof(stage.exp, dtg.stage, sigs_dev) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_PFC_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)
  gof_res[[paste0(stage, "_gof_res")]] <- gof

  ## MultiBrain

  dtg.MB <- dtg_res$MB %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    dplyr::filter(Sample %in% colnames(stage.exp)) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-contains(".Sig"))

  gof.MB <- write.gof(stage.exp, dtg.MB, sigsBrain$MB) %>%
    as.data.frame() %>%
    mutate(sigs = paste0(stage, "_MB_sig")) %>%
    rownames_to_column("SampleID") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\X", replacement = "") %>%
    mutate_at(.vars = "SampleID", .funs = gsub, pattern = "\\.", replacement = "-") %>%
    left_join(stage.md)

  gof_res[[paste0(stage, "_gof_MB")]] <- gof.MB
}
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 11296 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in compute.nmae(a, b): Vector of true values contains 8572 NA !!! NA
#;-) excluded
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Warning in max(X, na.rm = T): no non-missing arguments to max; returning -Inf
#;-) Warning in min(X, na.rm = T): no non-missing arguments to min; returning Inf
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
#;-) Joining, by = "SampleID"
```

``` r
gof_res %>%
  do.call(rbind, .) %>%
  ggplot(aes(sigs, r, fill = sigs)) +
  geom_violin() +
  facet_wrap(~Dev.stage) +
  theme(axis.text.x = element_blank())
#;-) Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).
```

![](https://i.imgur.com/rBotrd1.png)<!-- -->

<details style="margin-bottom:10px;">
<summary>
Session info
</summary>

``` r
sessioninfo::session_info()
#;-) ─ Session info ───────────────────────────────────────────────────────────────
#;-)  setting  value
#;-)  version  R version 4.2.2 Patched (2022-11-10 r83330)
#;-)  os       Ubuntu 22.04.1 LTS
#;-)  system   x86_64, linux-gnu
#;-)  ui       X11
#;-)  language (EN)
#;-)  collate  en_AU.UTF-8
#;-)  ctype    en_AU.UTF-8
#;-)  tz       Australia/Adelaide
#;-)  date     2023-01-30
#;-)  pandoc   2.19.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
#;-) 
#;-) ─ Packages ───────────────────────────────────────────────────────────────────
#;-)  package              * version   date (UTC) lib source
#;-)  AnnotationDbi        * 1.58.0    2022-04-26 [1] Bioconductor
#;-)  assertthat             0.2.1     2019-03-21 [1] CRAN (R 4.0.3)
#;-)  backports              1.4.1     2021-12-13 [1] CRAN (R 4.2.0)
#;-)  Biobase              * 2.56.0    2022-04-26 [1] Bioconductor
#;-)  BiocFileCache          2.4.0     2022-04-26 [1] Bioconductor
#;-)  BiocGenerics         * 0.42.0    2022-04-26 [1] Bioconductor
#;-)  BiocIO                 1.6.0     2022-04-26 [1] Bioconductor
#;-)  BiocParallel           1.30.4    2022-10-11 [1] Bioconductor
#;-)  biomaRt              * 2.52.0    2022-04-26 [1] Bioconductor
#;-)  Biostrings             2.64.1    2022-08-18 [1] Bioconductor
#;-)  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.2.2)
#;-)  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.0.3)
#;-)  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.0.5)
#;-)  blob                   1.2.3     2022-04-10 [1] CRAN (R 4.2.0)
#;-)  broom                  1.0.2     2022-12-15 [1] CRAN (R 4.2.2)
#;-)  cachem                 1.0.6     2021-08-19 [1] CRAN (R 4.1.1)
#;-)  cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.0.3)
#;-)  class                  7.3-20    2022-01-13 [1] CRAN (R 4.2.0)
#;-)  cli                    3.6.0     2023-01-09 [1] CRAN (R 4.2.2)
#;-)  codetools              0.2-18    2020-11-04 [1] CRAN (R 4.0.3)
#;-)  colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.2.2)
#;-)  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.2.1)
#;-)  curl                   4.3.3     2022-10-06 [1] CRAN (R 4.2.2)
#;-)  data.table           * 1.14.6    2022-11-16 [1] CRAN (R 4.2.2)
#;-)  DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.2.0)
#;-)  dbplyr                 2.2.1     2022-06-27 [1] CRAN (R 4.2.0)
#;-)  DeconRNASeq          * 1.38.0    2022-04-26 [1] Bioconductor
#;-)  DelayedArray           0.22.0    2022-04-26 [1] Bioconductor
#;-)  digest                 0.6.31    2022-12-11 [1] CRAN (R 4.2.2)
#;-)  dplyr                * 1.0.10    2022-09-01 [1] CRAN (R 4.2.1)
#;-)  dtangle              * 2.0.9     2019-12-01 [1] CRAN (R 4.2.2)
#;-)  dtw                    1.23-1    2022-09-19 [1] CRAN (R 4.2.2)
#;-)  DTWBI                * 1.1       2018-07-11 [1] CRAN (R 4.2.2)
#;-)  e1071                  1.7-12    2022-10-24 [1] CRAN (R 4.2.2)
#;-)  edgeR                * 3.38.4    2022-08-07 [1] Bioconductor
#;-)  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.1.0)
#;-)  entropy                1.3.1     2021-10-02 [1] CRAN (R 4.2.2)
#;-)  evaluate               0.19      2022-12-13 [1] CRAN (R 4.2.2)
#;-)  fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.2.2)
#;-)  farver                 2.1.1     2022-07-06 [1] CRAN (R 4.2.1)
#;-)  fastmap                1.1.0     2021-01-25 [1] CRAN (R 4.0.3)
#;-)  filelock               1.0.2     2018-10-05 [1] CRAN (R 4.1.0)
#;-)  forcats              * 0.5.2     2022-08-19 [1] CRAN (R 4.2.1)
#;-)  fs                     1.6.0     2023-01-23 [1] CRAN (R 4.2.2)
#;-)  gargle                 1.2.1     2022-09-08 [1] CRAN (R 4.2.1)
#;-)  generics               0.1.3     2022-07-05 [1] CRAN (R 4.2.1)
#;-)  GenomeInfoDb         * 1.32.4    2022-09-06 [1] Bioconductor
#;-)  GenomeInfoDbData       1.2.8     2022-06-30 [1] Bioconductor
#;-)  GenomicAlignments      1.32.1    2022-07-24 [1] Bioconductor
#;-)  GenomicFeatures      * 1.48.4    2022-09-20 [1] Bioconductor
#;-)  GenomicRanges        * 1.48.0    2022-04-26 [1] Bioconductor
#;-)  ggplot2              * 3.4.0     2022-11-04 [1] CRAN (R 4.2.2)
#;-)  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.2.0)
#;-)  googledrive            2.0.0     2021-07-08 [1] CRAN (R 4.1.0)
#;-)  googlesheets4          1.0.1     2022-08-13 [1] CRAN (R 4.2.1)
#;-)  gtable                 0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
#;-)  haven                  2.5.1     2022-08-22 [1] CRAN (R 4.2.1)
#;-)  highr                  0.10      2022-12-22 [1] CRAN (R 4.2.2)
#;-)  hms                    1.1.2     2022-08-19 [1] CRAN (R 4.2.1)
#;-)  htmltools              0.5.4     2022-12-07 [1] CRAN (R 4.2.2)
#;-)  httr                   1.4.4     2022-08-17 [1] CRAN (R 4.2.1)
#;-)  IRanges              * 2.30.1    2022-08-18 [1] Bioconductor
#;-)  jsonlite               1.8.4     2022-12-06 [1] CRAN (R 4.2.2)
#;-)  KEGGREST               1.36.3    2022-07-12 [1] Bioconductor
#;-)  knitr                  1.41      2022-11-18 [1] CRAN (R 4.2.2)
#;-)  labeling               0.4.2     2020-10-20 [1] CRAN (R 4.0.3)
#;-)  lattice                0.20-45   2021-09-22 [1] CRAN (R 4.1.1)
#;-)  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.2.2)
#;-)  limma                * 3.52.4    2022-09-27 [1] Bioconductor
#;-)  limSolve             * 1.5.6     2019-11-14 [1] CRAN (R 4.2.2)
#;-)  locfit                 1.5-9.7   2023-01-02 [1] CRAN (R 4.2.2)
#;-)  lpSolve                5.6.17    2022-10-10 [1] CRAN (R 4.2.2)
#;-)  lsa                    0.73.3    2022-05-09 [1] CRAN (R 4.2.2)
#;-)  lubridate              1.9.0     2022-11-06 [1] CRAN (R 4.2.2)
#;-)  magrittr             * 2.0.3     2022-03-30 [1] CRAN (R 4.2.0)
#;-)  MASS                   7.3-58.1  2022-08-03 [1] CRAN (R 4.2.1)
#;-)  Matrix                 1.5-3     2022-11-11 [1] CRAN (R 4.2.2)
#;-)  MatrixGenerics       * 1.8.1     2022-06-26 [1] Bioconductor
#;-)  matrixStats          * 0.63.0    2022-11-18 [1] CRAN (R 4.2.2)
#;-)  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.2.0)
#;-)  mime                   0.12      2021-09-28 [1] CRAN (R 4.1.1)
#;-)  ModelMetrics         * 1.2.2.2   2020-03-17 [1] CRAN (R 4.2.2)
#;-)  modelr                 0.1.10    2022-11-11 [1] CRAN (R 4.2.2)
#;-)  munsell                0.5.0     2018-06-12 [1] CRAN (R 4.0.3)
#;-)  pander               * 0.6.5     2022-03-18 [1] CRAN (R 4.2.0)
#;-)  pcaMethods           * 1.88.0    2022-04-26 [1] Bioconductor
#;-)  pillar                 1.8.1     2022-08-19 [1] CRAN (R 4.2.1)
#;-)  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.0.3)
#;-)  plyr                 * 1.8.8     2022-11-11 [1] CRAN (R 4.2.2)
#;-)  png                    0.1-8     2022-11-29 [1] CRAN (R 4.2.2)
#;-)  preprocessCore       * 1.58.0    2022-04-26 [1] Bioconductor
#;-)  prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.0.3)
#;-)  progress               1.2.2     2019-05-16 [1] CRAN (R 4.0.3)
#;-)  proxy                  0.4-27    2022-06-09 [1] CRAN (R 4.2.0)
#;-)  purrr                * 1.0.1     2023-01-10 [1] CRAN (R 4.2.2)
#;-)  quadprog               1.5-8     2019-11-20 [1] CRAN (R 4.1.0)
#;-)  R.cache                0.16.0    2022-07-21 [1] CRAN (R 4.2.1)
#;-)  R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
#;-)  R.oo                   1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
#;-)  R.utils                2.12.2    2022-11-11 [1] CRAN (R 4.2.2)
#;-)  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.1.1)
#;-)  rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.0.5)
#;-)  Rcpp                   1.0.10    2023-01-22 [1] CRAN (R 4.2.2)
#;-)  RCurl                  1.98-1.9  2022-10-03 [1] CRAN (R 4.2.1)
#;-)  readr                * 2.1.3     2022-10-01 [1] CRAN (R 4.2.1)
#;-)  readxl                 1.4.1     2022-08-17 [1] CRAN (R 4.2.1)
#;-)  reprex                 2.0.2     2022-08-17 [1] CRAN (R 4.2.1)
#;-)  reshape2               1.4.4     2020-04-09 [1] CRAN (R 4.0.3)
#;-)  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.2.0)
#;-)  rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.2.0)
#;-)  rlang                  1.0.6     2022-09-24 [1] CRAN (R 4.2.1)
#;-)  rlist                  0.4.6.2   2021-09-03 [1] CRAN (R 4.2.2)
#;-)  rmarkdown              2.19      2022-12-15 [1] CRAN (R 4.2.2)
#;-)  Rsamtools              2.12.0    2022-04-26 [1] Bioconductor
#;-)  RSQLite                2.2.20    2022-12-22 [1] CRAN (R 4.2.2)
#;-)  rstudioapi             0.14      2022-08-22 [1] CRAN (R 4.2.1)
#;-)  rtracklayer            1.56.1    2022-06-23 [1] Bioconductor
#;-)  rvest                  1.0.3     2022-08-19 [1] CRAN (R 4.2.1)
#;-)  S4Vectors            * 0.34.0    2022-04-26 [1] Bioconductor
#;-)  scales                 1.2.1     2022-08-20 [1] CRAN (R 4.2.1)
#;-)  sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.2.0)
#;-)  SnowballC              0.7.0     2020-04-01 [1] CRAN (R 4.2.2)
#;-)  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.2.2)
#;-)  stringr              * 1.5.0     2022-12-02 [1] CRAN (R 4.2.2)
#;-)  styler                 1.8.1     2022-11-07 [1] CRAN (R 4.2.2)
#;-)  SummarizedExperiment * 1.26.1    2022-04-29 [1] Bioconductor
#;-)  tibble               * 3.1.8     2022-07-22 [1] CRAN (R 4.2.1)
#;-)  tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.2.2)
#;-)  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.2.2)
#;-)  tidyverse            * 1.3.2     2022-07-18 [1] CRAN (R 4.2.1)
#;-)  timechange             0.1.1     2022-11-04 [1] CRAN (R 4.2.2)
#;-)  tzdb                   0.3.0     2022-03-28 [1] CRAN (R 4.2.0)
#;-)  utf8                   1.2.2     2021-07-24 [1] CRAN (R 4.1.0)
#;-)  vctrs                  0.5.2     2023-01-23 [1] CRAN (R 4.2.2)
#;-)  withr                  2.5.0     2022-03-03 [1] CRAN (R 4.2.0)
#;-)  xfun                   0.36      2022-12-21 [1] CRAN (R 4.2.2)
#;-)  XML                    3.99-0.13 2022-12-04 [1] CRAN (R 4.2.2)
#;-)  xml2                   1.3.3     2021-11-30 [1] CRAN (R 4.2.0)
#;-)  XVector                0.36.0    2022-04-26 [1] Bioconductor
#;-)  yaml                   2.3.6     2022-10-18 [1] CRAN (R 4.2.2)
#;-)  zlibbioc               1.42.0    2022-04-26 [1] Bioconductor
#;-) 
#;-)  [1] /usr/local/lib/R/site-library
#;-)  [2] /usr/lib/R/site-library
#;-)  [3] /usr/lib/R/library
#;-) 
#;-) ──────────────────────────────────────────────────────────────────────────────
```

</details>
