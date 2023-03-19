###
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "viridis", 
         "sunburstR", "pheatmap")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})


## Functions 
add_feature = function(feature_column, features){
    as.vector(sapply(feature_column, function(x){
        names(features)[sapply(features, function(f) x %in% f)]})) 
}


num_to_round = function(age){
    if (is.na(age)) {
        NaN
    } else if (age >= 2) {
        paste0(round(age), " yrs")
    } else if (age < 0) {
        paste0(round(age * 52 + 40), " pcw")
    } else if (age >= 0 & age < 2) {
        paste0(round(age * 12), " mos")
    }
}



## This calculates number of samples per age interval and number of samples 
summarise_stats = function(x, dataset)
{
    age = table(x$AgeInterval) %>% melt()
    age$Type = c(paste(dataset, "Sample", sep = "_"))
    individuals = x %>% group_by(AgeInterval, DonorID) %>% 
        summarise(n = n()) %>% 
        as.data.frame() 
    individuals = table(individuals$AgeInterval) %>% melt()
    individuals$Type = c(paste("Individual", dataset, sep = "_"))
    
    age = rbind(age, individuals)
    colnames(age) = c("AgeInterval", "n", "Type")
    age$dataset = c(as.character(dataset))
    return(age)
    
}



### Heatmap showing which samples have what characteristics 

of_interest = c("SampleID", "mito_Rate", "rRNA_rate", "RIN", "AgeNumeric", "Sex", 
                "AgeInterval", "Regions", "StructureAcronym", "Ethnicity", "Diagnosis", 
                "CellType", "PMI", "DonorID", "TotalNReads")




#### annotation data.

bspan.annot = read.csv("../annotations/BrainSpan-metadata-annot.csv", header = TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"BrainSpan" =`Include..Yes.No....Interest`) 


bseq.annot = read.csv("../annotations/BrainSeq-metadata-annot.csv", header = TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"BrainSeq" =`Include..Yes.No....Interest`) 

gtex.annot = read.csv("../annotations/GTEx-metadata-annot.csv", header = TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"GTEx" =`Include..Yes.No....Interest`) 

pe.annot = read.csv("../annotations/PsychEncode-metadata-annot.csv", header = TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"PsychEncode" =`Include..Yes.No....Interest`) 

vel.annot = read.csv("../annotations/Velmeshev-annot.csv", header=TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"Velmeshev" =`Include..Yes.No....Interest`) 

hca.annot = read.csv("../annotations/HCA-annot.csv", header=TRUE) %>% 
    dplyr::filter(`Include..Yes.No....Interest` == "Yes") %>%
    dplyr::filter(BITColumnName %in% of_interest) %>% 
    dplyr::select(BITColumnName,"HCA" =`Include..Yes.No....Interest`) 


dataset_annot = merge(gtex.annot, pe.annot,  by = "BITColumnName", all = TRUE)
dataset_annot = merge(bspan.annot, dataset_annot,  by = "BITColumnName", all = TRUE)
dataset_annot = merge(bseq.annot, dataset_annot,  by = "BITColumnName", all = TRUE)
dataset_annot = merge(hca.annot, dataset_annot,  by = "BITColumnName", all = TRUE)
dataset_annot = merge(vel.annot, dataset_annot,  by = "BITColumnName", all = TRUE)

dataset_annot= dataset_annot %>% 
    mutate_all(~replace(., is.na(.), 0)) 

dataset_annot = dataset_annot[-13,]

dataset_annot %<>% rownames_to_column("Random") %>%
    column_to_rownames("BITColumnName") 

dataset_annot %<>% dplyr::select(-Random)
names = rownames(dataset_annot)
names
dataset_annot = data.frame(apply(dataset_annot,2, function(x) as.numeric(gsub("Yes", 1,x))))
rownames(dataset_annot) = names

dataset_annot$labels = c("Phenotype", "Phenotype", "Sample", "Phenotype", 
                         "Phenotype", "Phenotype", "Seq", "Seq", "Sample",
                         "Seq", "Seq", "Sample", "Phenotype", "Sample", "Seq")

dataset_annot %<>% 
    arrange(labels)

dataset_annot
md = dataset_annot[,-7] %>% t() %>%
    pheatmap(cellheight = 15, 
             cellwidth = 15, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE, 
             border_color = "#FFF8F9", 
             color = c("#FFF8F9", "black"))

ggsave(file = "../../Results/exploratory/metadata.svg", md, 
       height = 7.5928, width = 19.9329, units = "cm")
#### Making an age interval and sample number graphs 

age_interval_stats = list()
directory = file.path("/home/neuro/Documents/Brain_integrative_transcriptome/Results/Formatted")
pattern = "/home/neuro/Documents/Brain_integrative_transcriptome/Results/Formatted/"

for (f in directory){
    md = list.files(f, full.names = TRUE, pattern = "\\-metadata.csv$")
    
    for (j in md){
        ct_file = read.csv(j, header= TRUE)
        dataset = gsub(pattern, "", j)
        dataset = gsub("\\-metadata.csv","", dataset)
        message("Now calculating statistics for ", dataset)
        stats = summarise_stats(ct_file, dataset)
        age_interval_stats[[paste(dataset)]] = stats
    }
    
}




bulk_plot = age_interval_stats %>% 
    do.call(rbind, .) %>% 
    mutate(AgeInterval = factor(AgeInterval, levels = c("4-7pcw", "8-9pcw",
                                                        "10-12pcw", "13-15pcw", "16-18pcw",
                                                        "19-24pcw", "25-38pcw", "0-5mos",
                                                        "6-18mos", "19mos-5yrs", "6-11yrs",
                                                        "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs",
                                                        "50-59yrs", "60-69yrs", "70-79yrs", "80-89yrs", "90-99yrs")), 
           Type = factor(Type, levels = c("BrainSeq_Sample", "Individual_BrainSeq", 
                                                      "BrainSpan_Sample", "Individual_BrainSpan", 
                                                      "GTEx_Sample", "Individual_GTEx", 
                                                      "PsychEncode_Sample", "Individual_PsychEncode"))) %>%
    drop_na() %>%
    ggplot(aes(x= AgeInterval, y = n, fill =Type)) +
        geom_bar(position = "dodge",stat= "identity") + 
        facet_grid(dataset~AgeInterval, scales = "free") + xlab("") + ylab("")  + theme_bw() +
        theme(legend.position = "none", axis.text.x=element_blank(),
              strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
              panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) + 
        theme(strip.text = element_text(colour = 'white')) +
        scale_fill_manual(values = c("#F75E5E", "#FFC6BD",
                                     "#49165E", "#EBBAFF",
                                     "#78A2EB", "#36466A",
                                     "#F9AD79", "#FF5F0F")) + 
        theme(
            panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
            plot.background = element_rect(fill = "transparent",colour =NA)
        )


bulk_plot

500ggsave(file = "../../Results/exploratory/bulk_dist_poster.svg", bulk_plot, 
       height = 14.9624, width = 39.3501, units = "cm")
    
### total number of donors and samples 
    
## Max value is 2642

total_samples_bulk = age_interval_stats %>% 
    do.call(rbind, .) %>%
    drop_na() %>%
    group_by(Type) %>% 
    summarise(sum = sum(n, na.rm = TRUE)) %>%
    as.data.frame()

brainseq_total = total_samples_bulk %>%
    drop_na() %>%
    dplyr::filter(grepl("BrainSeq", Type)) %>% 
    mutate(Type = factor(Type, levels = c("Individual_BrainSeq","BrainSeq_Sample", 
                                    "Individual_BrainSpan", "BrainSpan_Sample",
                                    "Individual_GTEx", "GTEx_Sample",
                                    "Individual_PsychEncode", "PsychEncode_Sample"))) %>%
    ggplot(aes(x= Type, y = sum, fill =Type)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c( "#FFC6BD","#F75E5E",
                                 "#49165E", "#EBBAFF",
                                 "#78A2EB", "#36466A",
                                 "#F9AD79", "#FF5F0F")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ylim(0, 2700) + coord_flip()

brainseq_total

ggsave(file = "../../Results/exploratory/BrainSeq_total.svg", brainseq_total, width = 1.98, 
       height = 1.89,, units = "in")
   

brainSpan_total = total_samples_bulk %>%
    drop_na() %>%
    dplyr::filter(grepl("BrainSpan", Type)) %>% 
    mutate(Type = factor(Type, levels = c("Individual_BrainSeq","BrainSeq_Sample", 
                                          "Individual_BrainSpan", "BrainSpan_Sample",
                                          "Individual_GTEx", "GTEx_Sample",
                                          "Individual_PsychEncode", "PsychEncode_Sample"))) %>%
    ggplot(aes(x= Type, y = sum, fill =Type)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c("#EBBAFF", "#49165E")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ylim(0, 2700) + coord_flip()


ggsave(file = "../../Results/exploratory/BrainSpan_total.svg", brainSpan_total, width = 1.98, 
       height = 1.89, units = "in")


### GTex
gtex_total = total_samples_bulk %>%
    drop_na() %>%
    dplyr::filter(grepl("GTEx", Type)) %>% 
    mutate(Type = factor(Type, levels = c("Individual_BrainSeq","BrainSeq_Sample", 
                                          "Individual_BrainSpan", "BrainSpan_Sample",
                                          "Individual_GTEx", "GTEx_Sample",
                                          "Individual_PsychEncode", "PsychEncode_Sample"))) %>%
    ggplot(aes(x= Type, y = sum, fill =Type)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c( "#36466A", "#78A2EB")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ylim(0, 2700) + coord_flip()


ggsave(file = "../../Results/exploratory/GTEx_total.svg",gtex_total,  width = 1.98, 
       height = 1.89, units = "in")


### PsychEncode 

pe_total = total_samples_bulk %>%
    drop_na() %>%
    dplyr::filter(grepl("Psych", Type)) %>% 
    mutate(Type = factor(Type, levels = c("Individual_BrainSeq","BrainSeq_Sample", 
                                          "Individual_BrainSpan", "BrainSpan_Sample",
                                          "Individual_GTEx", "GTEx_Sample",
                                          "Individual_PsychEncode", "PsychEncode_Sample"))) %>%
    ggplot(aes(x= Type, y = sum, fill =Type)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c( "#FF5F0F",  "#F9AD79")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ylim(0, 2700) + coord_flip()


pe_total
ggsave(file = "../../Results/exploratory/PE_total.svg",pe_total, width = 1.98, 
       height = 1.89, units = "in")



## Single-cell age distribution


hca = read.csv("../../../BrainData/single-cell/hca/HCA-metadata_fixed.csv", header=TRUE)
vel = read.csv("../../../BrainData/single-cell/velmeshev/Velmeshev-metadata.csv", header= TRUE)

dim(vel)
hca %<>% mutate(Age_rounded = as.character(sapply(na.omit(.$AgeNumeric), num_to_round))) %>%
    mutate(AgeInterval = add_feature(.$Age_rounded, age_intervals))

vel %<>% 
    mutate(Age_rounded = as.character(sapply(na.omit(.$AgeNumeric), num_to_round))) %>%
    mutate(AgeInterval = add_feature(.$Age_rounded, age_intervals))


vel_age = table(vel$AgeInterval) %>% melt()
vel_age$dataset = "Velmeshev"

vel %>% group_by(AgeInterval,  MajorCelltype) %>% 
    summarise(n = n()) %>% as.data.frame

hca_age = table(hca$AgeInterval) %>% melt()
hca_age$dataset = "HCA"
sc_age_dist = rbind(hca_age,vel_age) %>% 
    add_row(Var1 = c("4-7pcw", "8-9pcw", "10-12pcw", "13-15pcw", "16-18pcw", "19-24pcw", 
                     "25-38pcw", "0-5mos","6-18mos",  "30-39yrs", "60-69yrs",
                     "70-79yrs", "80-89yrs", "90-99yrs"), value = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
            dataset = c("Velmeshev")) %>%
    mutate(Var1 = factor(Var1, levels = c("4-7pcw", "8-9pcw", "10-12pcw", "13-15pcw", "16-18pcw", "19-24pcw", 
                                          "25-38pcw", "0-5mos","6-18mos", "19mos-5yrs", "6-11yrs","12-19yrs", 
                                          "20-29yrs", "30-39yrs", "40-49yrs","50-59yrs", "60-69yrs",
                                          "70-79yrs", "80-89yrs", "90-99yrs"))) %>%
    drop_na() %>%
    ggplot(aes(x= Var1, y = value, fill =dataset)) +
    geom_bar(position = "dodge",stat= "identity", width = 0.5) + 
    facet_grid(dataset~Var1, scales = "free") + xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c("#07473C","#B82828")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA)
    )

sc_age_dist

ggsave(file = "../../Results/exploratory/sc_dist_poster.svg", sc_age_dist, 
       height = 8.0819, width = 39.3501, units = "cm")


## Number of total cells 

vel_total= data.frame(dataset = c("Vel", "HCA"), 
           val = c(nrow(vel), nrow(hca))) %>%
    dplyr::filter(grepl("Vel", dataset)) %>% 
    ggplot(aes(x= dataset, y = val, fill =dataset)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c( "#B82828")) + 
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour =NA), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + ylim(0, 85000) + coord_flip()

ggsave(file = "../../Results/exploratory/vel_total.svg",vel_total, width = 1.98, 
       height = 1.89, units = "in")


data.frame(dataset = c("Vel", "HCA"), 
           val = c(nrow(vel), nrow(hca))) %>%
    ggplot(aes(x= dataset, y = val, fill =dataset)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw()  + ylim(0, 85000) + coord_flip()


hca_total = data.frame(dataset = c("Vel", "HCA"), 
                      val = c(nrow(vel), nrow(hca))) %>%
    dplyr::filter(grepl("HCA", dataset)) %>% 
    ggplot(aes(x= dataset, y = val, fill =dataset)) +
    geom_bar(position = "dodge",stat= "identity", width = .5) +  xlab("") + ylab("")  + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(),
          strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
          panel.border = element_blank()) + 
    theme(strip.text = element_text(colour = 'white')) +
    scale_fill_manual(values = c("#07473C"))  + ylim(0, 85000) + coord_flip()

hca_total
ggsave(file = "../../Results/exploratory/hca_total.svg",hca_total, width = 1.98, 
       height = 1.89, units = "in")
## 20000, 40,000, 60,000, 80,000


vel %<>% 
    mutate(StructureAcronym = gsub("PFC", "DLPFC", .$StructureAcronym)) %>% 
    mutate(Regions = add_feature(.$StructureAcronym, regions))

hca %<>% 
    mutate(Regions = add_feature(.$region, regions))


### Regions dataset 

#### snRNA-seq 
regions_data_sc= vel %>% dplyr::select(CellID, StructureAcronym, Regions)
regions_data_sc$dataset = "Velmeshev"

colnames(regions_data_sc)[1] = c("SampleID")
head(hca)
regions_data_sc = hca %>% 
    dplyr::select("SampleID" = "sample_name", "StructureAcronym" = "region", Regions) %>%
    mutate(dataset = c("HCA")) %>% 
    rbind(., regions_data_sc)
    
table(regions_data_sc$dataset)    


write.csv(regions_data_sc,"../../Results/exploratory/regions_sc.csv")
## Bulk data

## Bspan
regions_data

regions_data = read.csv("../../../BrainData/Bulk/BrainSpan/Formatted/BrainSpan-metadata.csv", header=TRUE) %>% 
    dplyr::select(SampleID, StructureAcronym, Regions) %>% 
    mutate(dataset = c("BrainSpan")) %>% 
    rbind(., regions_data)
    
regions_data = read.csv("../../Results/Formatted/BrainSeq-metadata.csv", header = TRUE) %>%
    dplyr::select(SampleID, StructureAcronym, Regions) %>% 
    mutate(dataset = c("BrainSeq")) %>% 
    rbind(., regions_data)


regions_data = read.csv("../../Results/Formatted/GTEx-metadata.csv", header = TRUE) %>%
    dplyr::select(SampleID, StructureAcronym, Regions) %>% 
    mutate(dataset = c("GTEx")) %>% 
    rbind(., regions_data)


regions_data = read.csv("../../Results/Formatted/PsychEncode-metadata.csv", header = TRUE) %>%
    dplyr::select(SampleID, StructureAcronym, Regions) %>% 
    mutate(StructureAcronym = gsub("PFC", "DLPFC", .$StructureAcronym)) %>%
    mutate(dataset = c("PsychEncode")) %>% 
    rbind(., regions_data)


table(regions_data$StructureAcronym)

write.csv(regions_data, "../../Results/exploratory/regions_bulk_sc.csv")



## snRNA-seq plots - circular barplots


### example for circular barplot 

hike_data <- readr::read_rds("../../Data/hike_data.rds")
hike_data

hike_data$region <- as.factor(word(hike_data$location, 1, sep = " -- "))
hike_data$length_num <- as.numeric(sapply(strsplit(hike_data$length, " "), "[[", 1))

hike_data


plot_df <- hike_data %>%
    group_by(region) %>%
    summarise(
        sum_length = sum(length_num),
        mean_gain = mean(as.numeric(gain)),
        n = n()
    ) %>%
    mutate(mean_gain = round(mean_gain, digits = 0))


plot_df


plt <- ggplot(plot_df) +
    # Make custom panel grid
    geom_hline(
        aes(yintercept = y), 
        data.frame(y = c(0:3) * 1000),
        color = "lightgrey"
    ) + 
    # Add bars to represent the cumulative track lengths
    # str_wrap(region, 5) wraps the text so each line has at most 5 characters
    # (but it doesn't break long words!)
    geom_col(
        aes(
            x = reorder(str_wrap(region, 5), sum_length),
            y = sum_length,
            fill = n
        ),
        position = "dodge2",
        show.legend = TRUE,
        alpha = .9
    ) +
    
    # Add dots to represent the mean gain
    geom_point(
        aes(
            x = reorder(str_wrap(region, 5),sum_length),
            y = mean_gain
        ),
        size = 3,
        color = "gray12"
    ) +
    
    # Lollipop shaft for mean gain per region
    geom_segment(
        aes(
            x = reorder(str_wrap(region, 5), sum_length),
            y = 0,
            xend = reorder(str_wrap(region, 5), sum_length),
            yend = 3000
        ),
        linetype = "dashed",
        color = "gray12"
    ) + 
    
    # Make it circular!
    coord_polar()

plt


### plotting for HCA data 


hca.CT = hca %>% group_by(subclass, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    mutate(percentage = (n/ sum(n)) * 100) 

hca.CT$MajorCellType[hca.CT$MajorCelltype =="Astrocytes" | hca.CT$MajorCelltype == "Microglia" | 
                         hca.CT$MajorCelltype == "Oligodendrocytes" | hca.CT$MajorCelltype == "OPCs"] <- c("Neuroglia")


hca.CT$MajorCellType[hca.CT$MajorCelltype =="Vasculature"] <- c("Vasculature")
hca.CT$MajorCellType[grepl("Neurons", hca.CT$MajorCelltype, ignore.case = TRUE)] <- c("Neurons")
hca.CT

plot_CT  = ggplot(hca.CT) +
    # Make custom panel grid
    geom_hline(
        aes(yintercept = y), 
        data.frame(y = c(0:5) * 10),
        color = "lightgrey"
    ) + 
    # Add bars to represent the cumulative track lengths
    # str_wrap(region, 5) wraps the text so each line has at most 5 characters
    # (but it doesn't break long words!)
    geom_col(
        aes(
            x = reorder(str_wrap(subclass, 5), n),
            y = log10(n),
            fill = MajorCellType
        ),
        position = "dodge2",
        show.legend = TRUE,
        alpha = .9
    ) + 
    # Add dots to represent the mean gain
    geom_point(
        aes(
            x = reorder(str_wrap(subclass, 5),n),
            y = log10(n)
        ),
        size = 3,
        color = "gray12"
    ) +
    
    # Lollipop shaft for mean gain per region
    geom_segment(
        aes(
            x = reorder(str_wrap(subclass, 5), n),
            y = 0,
            xend = reorder(str_wrap(subclass, 5), n),
            yend = 3
        ),
        linetype = "dashed",
        color = "gray12"
    ) + 

    
    # Make it circular!
    coord_polar()

plot_CT=  plot_CT +

    # New fill and legend title for number of tracks per region
    # Make the guide for the fill discrete
    theme(
        # Remove axis ticks and text
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        # Use gray text for the region names
        axis.text.x = element_text(color = "gray12", size = 12),
        # Move the legend to the bottom
        legend.position = "bottom",
    ) + 
     scale_y_continuous(
         limits = c(-2, 4.5),
         expand = c(0, 0),
         breaks = c(0, 1000, 2000, 3000)
     )  + 
     scale_fill_manual(
         values = c("#FFFDC2", "#D94A6F", "#68C9CF")
     ) 
plot_CT

ggsave(file = "../../Results/exploratory/hca_CT.svg", width=7.11, height =6.90, 
       units = "in")



vel.ct = vel %>%  group_by(CellType, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    arrange(MajorCelltype, n)

vel.ct$MajorCellType[vel.ct$MajorCelltype =="Astrocytes" | vel.ct$MajorCelltype == "Microglia" | 
                                vel.ct$MajorCelltype == "Oligodendrocytes" | vel.ct$MajorCelltype == "OPCs"] <- c("Neuroglia")


vel.ct$MajorCellType[vel.ct$MajorCelltype =="Endothelia"] <- c("Vasculature")
vel.ct$MajorCellType[grepl("Neurons", vel.ct$MajorCelltype, ignore.case = TRUE)] <- c("Neurons")


plot_CT_vel  = ggplot(vel.ct) +
    # Make custom panel grid
    geom_hline(
        aes(yintercept = y), 
        data.frame(y = c(0:5) * 10),
        color = "lightgrey"
    ) + 
    # Add bars to represent the cumulative track lengths
    # str_wrap(region, 5) wraps the text so each line has at most 5 characters
    # (but it doesn't break long words!)
    geom_col(
        aes(
            x = reorder(str_wrap(CellType, 5), n),
            y = log10(n),
            fill = MajorCellType
        ),
        position = "dodge2",
        show.legend = TRUE,
        alpha = .9
    ) + 
    # Add dots to represent the mean gain
    geom_point(
        aes(
            x = reorder(str_wrap(CellType, 5),n),
            y = log10(n)
        ),
        size = 3,
        color = "gray12"
    ) +
    
    # Lollipop shaft for mean gain per region
    geom_segment(
        aes(
            x = reorder(str_wrap(CellType, 5), n),
            y = 0,
            xend = reorder(str_wrap(CellType, 5), n),
            yend = 3
        ),
        linetype = "dashed",
        color = "gray12"
    ) + 
    
    
    # Make it circular!
    coord_polar()

#plot_CT_vel =  
plot_CT_vel = plot_CT_vel +
    
    # New fill and legend title for number of tracks per region
    # Make the guide for the fill discrete
    theme(
        # Remove axis ticks and text
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        # Use gray text for the region names
        axis.text.x = element_text(color = "gray12", size = 12),
        # Move the legend to the bottom
        legend.position = "bottom",
    ) + 
    scale_y_continuous(
        limits = c(-2, 4.5),
        expand = c(0, 0),
        breaks = c(0, 1000, 2000, 3000)
    )  + 
    scale_fill_manual(
        values = c("#FFFDC2", "#D94A6F", "#68C9CF")
    ) 


plot_CT_vel




ggsave(file = "../../Results/exploratory/vel_CT.svg",plot_CT_vel, width=7.11, height =6.90, 
       units = "in")
