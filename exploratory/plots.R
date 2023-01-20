## Figures 
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "viridis")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})



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


directory = file.path("/Users/urwah/Documents/PhD/Brain_transcriptome/Data/FormattedData/Metadata")
pattern = "/Users/urwah/Documents/PhD/Brain_transcriptome/Data/FormattedData/Metadata/"

for (f in directory){
  md = list.files(f, full.names = TRUE, pattern = "\\-metadata.csv$")
  
  for (j in md){
    ct_file = read.csv(j, header= TRUE)
    dataset = gsub(pattern, "", j)
    dataset = gsub("\\-metadata.csv","", dataset)
    message("Now calculating statistics for ", dataset)
    stats = summarise_stats(ct_file)
    stats$Dataset = c(as.character(dataset))
  }
  
}



bspan.md = read.csv("../../Data/FormattedData/BrainSpan-metadata.csv")
bseq.md = read.csv("../../Data/FormattedData/BrainSeq-metadata.csv")
gtex.md = read.csv("../../Data/FormattedData/GTEx-metadata.csv")
pe.md = read.csv("../../Data/FormattedData/PsychEncode-metadata.csv", header=TRUE)


age = table(pe.md$AgeInterval) %>% melt()
age$Type = c("Sample")
individual = table(pe.md$AgeInterval) %>% melt()
individual$Type = c("Individual")
age = rbind(age, individual)
age$dataset = c("PsychEncoode")
colnames(age) = c("AgeInterval", "n", "Type", "dataset")
age


age = summarise_stats(gtex.md, "GTEx")
age
age.bspan = summarise_stats(bspan.md, "BrainSpan")
age.bseq = summarise_stats(bseq.md, "BrainSeq")

all_2 = rbind(age, age.bspan, age.bseq)

all_2 = rbind(all_2, age)
table(all_2$dataset)


write.csv(all_2, "../../Results/Metadata/ageIntervals.csv")

unique(all_2$AgeInterval)
all_2$AgeInterval = factor(all_2$AgeInterval, levels = c("4-7pcw", "8-9pcw",
                                                              "10-12pcw", "13-15pcw", "16-18pcw",
                                                              "19-24pcw", "25-38pcw", "0-5mos",
                                                              "6-18mos", "19mos-5yrs", "6-11yrs",
                                                              "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs",
                                                              "50-59yrs", "60-69yrs", "70-79yrs", "80-89yrs", "90-99yrs"))


all_2$Type = factor(all_2$Type, levels = c("BrainSeq_Sample", "Individual_BrainSeq", 
                                           "BrainSpan_Sample", "Individual_BrainSpan", 
                                           "GTEx_Sample", "Individual_GTEx", 
                                           "Sample", "Individual"))



bspan.md %>% 
  group_by(DonorID) %>% 
  summarise(n = n())

length(unique(bspan.md$DonorID))
length(unique(bspan.md$SampleID))
total_samples = data.frame(
  Dataset = c("BrainSpan", "BrainSeq", "GTEx", "PsychEncode"), 
  Samples = c(length(unique(bspan.md$SampleID)),length(unique(bseq.md$RNum)), 
              length(unique(gtex.md$SampleID)), length(unique(pe.md$SampleID))), 
  Individuals = c(length(unique(bspan.md$DonorID)), length(unique(bseq.md$DonorID)), 
                  length(unique(gtex.md$DonorID)), length(unique(pe.md$SampleID)))
)

total_samples$Dataset = factor(total_samples$Dataset, 
                               levels = c("PsychEncode",
                                           "GTEx",
                                           "BrainSpan", 
                                           "BrainSeq"))


total_samples %<>% melt %>% mutate(fill_col = paste(Dataset, variable, sep = "_"))
write.csv(total_samples, "../../Results/Metadata/total_samples.csv")

total= total_samples %>% 
  melt() %>%
  mutate(fill_col = paste(Dataset, variable, sep = "_")) %>%
  ggplot(aes(x = Dataset, y =value, fill= fill_col)) +
  geom_bar(position = "dodge",stat= "identity", color = "#B69E96") + 
  xlab("") + ylab("") + theme_bw() +
  theme(legend.position = "none",
        strip.background =element_rect( color = "#E1DFDB"), 
        axis.text.x = element_text(face="bold", size = 10),
        axis.text.y = element_text(face="bold", size =10),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
  theme(strip.text = element_text(colour = 'white')) +
  scale_fill_manual(values = c("#92373F","#A5635E", 
                                        "#A76972","#BB9992", 
                                        "#E89787", "#CA756E", 
                                        "#AD9783", "#C1CCA6")) +
                                          theme(panel.background = 
                                                  element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                                plot.background = element_rect(fill = "transparent",colour = NA)) + coord_flip()
  

ggsave(file="../../Results/Metadata/TotalSamples.svg", plot = total, width=6.18, height = 5.6, 
       units= "in", device="svg")


## regional plots 

regions = read.csv("../../Results/Metadata/bulk_regions.csv", header=TRUE)

regions$AgeInterval =factor(regions$AgeInterval, levels = c("4-7pcw", "8-9pcw",
                                                                               "10-12pcw", "13-15pcw", "16-18pcw",
                                                                               "19-24pcw", "25-38pcw", "0-5mos",
                                                                               "6-18mos", "19mos-5yrs", "6-11yrs",
                                                                               "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs",
                                                                               "50-59yrs", "60-69yrs", "70-79yrs", "80-89yrs", "90-99yrs"))

regions %>% 
  drop_na() %>%
  group_by(Dataset,AgeInterval,Regions) %>% 
  summarise(n = n()) %>% as.data.frame() %>%
  ggplot(aes(x=AgeInterval, y=n, fill=Regions)) +
  geom_bar(stat= "identity", color = "black") + 
  facet_grid(Dataset ~ AgeInterval,scales = "free") + xlab("") + ylab("")  +
  theme_classic() + 
  theme(strip.text = element_text(colour = 'white')) +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        strip.background =element_rect(fill="#B69E96", color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
  scale_fill_manual(values = c("#B15759","#C6A897", 
                                        "#B3C48F","#E4A28A")) +
                                          theme(panel.background = 
                                                  element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                                plot.background = element_rect(fill = "transparent",colour = NA))


B15759

regions$StructureAcronym = gsub("^PFC", "DLPFC", regions$StructureAcronym)
write.csv(regions, "../../Results/Metadata/bulk_regions.csv")
table(regions$StructureAcronym)

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
all_2 %<>% drop_na()
all
table(pe.md$AgeInterval)
ageInterval =all_2 %>% as.data.frame() %>% 
  drop_na() %>% 
  ggplot(aes(x= AgeInterval, y = n, fill =Type)) +
  geom_bar(position = "dodge",stat= "identity", color = "#B69E96") + 
  facet_grid(dataset ~ AgeInterval,scales = "free") + xlab("") + ylab("")  + theme_classic() +
  theme(legend.position = "none", axis.text.x=element_blank(),
        strip.background =element_rect(fill="#B69E96", color = "#E1DFDB"), 
        panel.border = element_rect(color = "#E1DFDB", fill = NA, size = 1)) + 
  theme(strip.text = element_text(colour = 'white')) +
  scale_fill_manual(values = c("#92373F","#A5635E", 
                                "#A76972","#BB9992",
                               "#E89787", "#CA756E",
                               "#AD9783", "#C1CCA6")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

ageInterval
ggsave(file="../../Results/Metadata/AgeInterval.svg", plot = ageInterval, width=16, height = 5.6, 
       units= "in", device="svg")


all_2 %>% drop_na() %>% 
  ggplot(aes(x = dataset, y =n, fill= Type)) +
  geom_bar(position = "dodge",stat= "identity") + 
  xlab("") + ylab("") + theme_bw() +
  theme(legend.position = "none",
        strip.background =element_rect( color = "#E1DFDB"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
  theme(strip.text = element_text(colour = 'white')) +
  scale_fill_manual(values = c("#92373F","#A5635E", 
                                        "#A76972","#BB9992", 
                                        "#E89787", "#CA756E", 
                                        "#AD9783", "#C1CCA6")) +
  theme(panel.background = 
          element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour = NA))

all_2 %>% head()


?ggsave
"#C58C8B"

c("")
theme(strip.background =element_rect(fill="red"))

bg_col = c("#8a616e")

















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



