## Figures 
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "viridis")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

#### Functions ####

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


directory = file.path("/home/neuro/Documents/BrainData/Bulk")

pattern = "/home/neuro/Documents/BrainData/Bulk/"

age_int = NULL
for (f in directory){
  md = list.files(f, full.names = TRUE, pattern = "\\-metadata.csv$", recursive = TRUE)
  
  for (j in md){
    ct_file = read.csv(j, header= TRUE)
  
    data = gsub(paste0(pattern, ), "", j)
    data = gsub("\\-metadata.csv","", data)
    message("Now calculating statistics for ", data)
    stats = summarise_stats(ct_file, as.character(data))
    age_int = rbind(age_int, stats)
    
  }
  
}

stats

bspan.md = read.csv("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv/Formatted/BrainSpan-metadata.csv")
bseq.md = read.csv("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-metadata.csv")
gtex.md = read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-metadata.csv")
pe.md = read.csv("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted/PsychEncode-metadata.csv", header=TRUE)


unique(gtex.md$Structure)
age = table(pe.md$AgeInterval) %>% melt()
age$Type = c("Sample")
individual = table(pe.md$AgeInterval) %>% melt()
individual$Type = c("Individual")
age = rbind(age, individual)
age$dataset = c("PsychEncoode")
colnames(age) = c("AgeInterval", "n", "Type", "dataset")

age_gtex = summarise_stats(gtex.md, "GTEx")

age.bspan = summarise_stats(bspan.md, "BrainSpan")
age.bseq = summarise_stats(bseq.md, "BrainSeq")

all = rbind(age, age.bspan, age.bseq, age_gtex)


unique(all$AgeInterval)
all$AgeInterval = factor(all$AgeInterval, levels = c("4-7pcw", "8-9pcw",
                                                              "10-12pcw", "13-15pcw", "16-18pcw",
                                                              "19-24pcw", "25-38pcw","39-40pcw" ,"0-5mos",
                                                              "6-18mos", "19mos-5yrs", "6-11yrs",
                                                              "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs",
                                                              "50-59yrs", "60-69yrs", "70-79yrs", "80-89yrs", "90-99yrs"))


all$Type = factor(all$Type, levels = c("BrainSeq_Sample", "Individual_BrainSeq", 
                                           "BrainSpan_Sample", "Individual_BrainSpan", 
                                           "GTEx_Sample", "Individual_GTEx", 
                                           "Sample", "Individual"))











## Single cell 

hca= read.csv("../../Data/FormattedData/Metadata/HCA-metadata.csv", header=TRUE)
vel = read.csv("../../Data/FormattedData/Metadata/Velmeshev-metadata.csv", header= TRUE)

hca
head(hca)
table(hca$AgeNumeric)
table(hca$MajorCelltype)
table(vel$MajorCelltype)
dim(hca) 
unique(all_2$Type)
dim(vel)
max(all$n)
all %<>% drop_na()
all %>% as.data.frame() %>% 
  ggplot(aes(x= AgeInterval, y = n, fill =Type)) +
  geom_bar(position = "dodge",stat= "identity") + 
  facet_grid(dataset ~ AgeInterval,scales = "free") + xlab("") + ylab("")  + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank(),
        strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
        panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) + 
  theme(strip.text = element_text(colour = 'white')) +
  scale_fill_manual(values = c("#573D5E","#88738D", 
                               "#795761", "#D2A8B6", 
                               "#7C5D4C", "#D8B3A2", 
                               "#485B59", "#97BABB")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

ggsave(file="../../Results/exploratory/bulk_distribution.svg", height = 6.48, width = 16.65, units = "in")






theme(strip.background =element_rect(fill="red"))

bg_col = c("#8a616e")









## Regional composition per age interval 

##### REGIONS 

region_dist = NULL 

datasets = c("bspan.md", "pe.md", "gtex.md", "bseq.md")

for (data in datasets){
  data = noquote(data)
  regions = table(data$Regions, data$AgeInterval) %>% melt
  colnames(regions) = c("Region", "AgeInterval", as.character(data))
  region_dist = full_join(regions_dist, regions)
}



bspan_regions <- table(bspan.md$Regions, bspan.md$AgeInterval) %>% melt()
colnames(bspan_regions) <- c("Region", "AgeInterval", "BrainSpan")
bspan_regions



bseq_regions <- table(bseq.md$Regions, bseq.md$AgeInterval) %>% melt()
colnames(bseq_regions) <- c("Region", "AgeInterval", "BrainSeq")


regions <- full_join(bspan_regions, bseq_regions)
regions
gtex_regions <- table(gtex.md$Regions, gtex.md$AgeInterval) %>% melt()
gtex_regions

colnames(gtex_regions) <- c("Region", "AgeInterval", "GTEx")
gtex_regions
regions <- full_join(regions, gtex_regions)


pe_regions <- table(pe.md$Regions, pe.md$AgeInterval) %>% melt()
colnames(pe_regions) <- c("Region", "AgeInterval", "PsychEncode")

regions <- full_join(regions, pe_regions)

regions[is.na(regions)] <- 0
regions %<>% melt()


regions

table(bspan.md$StructureAcronym) 

regions$AgeInterval <- factor(regions$AgeInterval, levels = c("4-7pcw", "8-9pcw", 
                                                              "10-12pcw", "13-15pcw", "16-18pcw", 
                                                              "19-24pcw", "25-38pcw","39-40pcw",  "0-5mos", 
                                                              "6-18mos", "19mos-5yrs", "6-11yrs", 
                                                              "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs", 
                                                              "50-59yrs", "60-69yrs", "70-79yrs", "80-89yrs", "90-99yrs"))




regions

regions %>% 
  summarise(sum = n(value[]))

age_all <- age %>%
  ggplot(aes(x=dataset, y= value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black", width=7, position = "dodge") + ylab("") +
  xlab("")  + theme_bw() +   theme(axis.text.x = element_blank(), axis.title.y = element_text(face="bold"), 
                                   legend.position = "none") + 
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "magma") + facet_grid(dataset ~ AgeInterval, scales = "free")


age_count <- age %>%
  ggplot(aes(x=fct_rev(dataset), y= value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + ylab("") +
  xlab("")  + theme_bw() +   theme(axis.text.x = element_text(face="bold", angle = 90)) + 
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "magma") + coord_flip() 


age_count
regions 

bs <- regions %>% dplyr::filter(variable == "BrainSpan")
table(bs$Region)
table(bs$AgeInterval)

regions

pdf(file= "BITHub_samp_regions.pdf", height=8, width = 16)

regions %>% 
  ggplot(aes(x = variable, y= value, fill = Region)) + 
  geom_bar(stat = "identity", position = "fill", width = 7) + ylab("") +
  xlab("")  + theme_bw() +   theme(axis.text.x = element_blank(), axis.title.y = element_text(face="bold")) +
   facet_grid(variable ~ AgeInterval, scales = "free") + 
  scale_fill_manual(values = c("#573D5E",  "#795761", "#7C5D4C", "#485B59")) +
  theme(legend.position = "none", axis.text.x=element_blank(),
        strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
        panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) +
  theme(strip.text = element_text(colour = 'white'))

region_plots

regions %>% drop_na() %>% as.data.frame() %>% 
  ggplot(aes(x = variable, y= value, fill = Region)) +
  geom_bar(position = "fill",stat= "identity", width =7, color = "black") + 
  facet_grid(variable ~ AgeInterval, scales = "free") + xlab("") + ylab("")  + theme_bw() +
  theme(axis.text.x=element_blank(), legend.position = "top",
        strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
        panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) + 
  theme(strip.text = element_text(colour = 'white')) +
 # scale_fill_manual(values = c("#795761","#573D5E",  "#7C5D4C", "#485B59")) +
 scale_fill_manual(values = c("#88738D", "#D2A8B6","#D8B3A2", "#97BABB")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

ggsave(file="../../Results/exploratory/bulk_region_composition-ex2.svg", height = 6.48, width = 16.65, units = "in")

## individuals and samples per region 


data.frame("Dataset" = c("BrainSpan", "BrainSeq", "GTEx", "PsychEncode"), 
           "Sample" = c(length(bspan.md$SampleID), length(bseq.md$SAMPLE_ID), length(gtex.md$SampleID),
                        length(pe.md$SampleID)), 
           "Individual"= c(length(unique(bspan.md$DonorID)), 
                           length(unique(bseq.md$DonorID)), 
                           length(unique(gtex.md$DonorID)), 
                           length(unique(pe.md$DonorID)))) %>% melt() %>% 
  mutate(fill_col = paste(.$Dataset, .$variable, sep = "_")) %>% 
  ggplot(aes(x= Dataset, y = value, fill =fill_col)) +
  geom_bar(position = "dodge",stat= "identity", color = "#E1DFDB") + xlab("") + ylab("")  + theme_bw() +
  theme(legend.position = "none",  strip.background =element_rect(fill="#AA9A9C", color = "#E1DFDB"), 
        panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1), 
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 10)) +
  scale_fill_manual(values = c("#573D5E","#88738D", 
                               "#795761", "#D2A8B6", 
                               "#7C5D4C", "#D8B3A2", 
                               "#485B59", "#97BABB")) + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(file="../../Results/exploratory/bulk_total.svg", height = 7.98, width = 6.27, units = "in")



length(unique(bseq.md$DonorID))
unique(bseq.md$DonorID)
unique(bseq.md$SAMPLE_ID)
length()

### scatter plot 

bh_results = read.csv("../../Results/chart.csv", header=TRUE)
min(bh_results$BrainSeq)


bh_results %>% ggplot(aes(x=BrainSpan, y =BrainSeq, color = highlighted)) + geom_point(size = 1, alpha = 0.5) + 
  theme_classic() + theme(legend.position = "none") + scale_color_manual(values = c("#D2A8B6", "#573D5E" )) +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) + ggtitle("Z-Score Transformed Mean Log2(Expression)") + 
  theme(plot.title = element_text(hjust = 0.5)) 


## variance partition 

bseq_varPart = read.csv("../../Data/FormattedData/BrainSpan/BrainSpan-varPart.csv", header= TRUE, row.names =1)
bseq_varPart["ENSG00000101040",] %>% melt() %>% 
  ggplot(aes(x= reorder(variable, -value), y = value)) + geom_bar(stat = "identity", fill = "#583f60") + 
  theme_bw() + ylab("Proportion of variance explained") + xlab("") + 
  theme(axis.text.y = element_text(hjust=1)) +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA)
  )



bseq_varPart["ENSG00000116273",]

dir = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")
outdir = file.path(dir, "Formatted")



num_to_round(-0.019178)
