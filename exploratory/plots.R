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



bspan.md = read.csv("../../../Data/FormattedData/BrainSpan/BrainSpan-metadata.csv")
bseq.md = read.csv("../../Data/FormattedData/BrainSeq/BrainSeq-metadata.csv")
gtex.md = read.csv("../../Data/FormattedData/Gtex/GTEx-metadata.csv")
pe.md = read.csv("../../Data/FormattedData/PsychEncode/PsychEncode-metadata.csv", header=TRUE)


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

unique(all$AgeInterval)
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
all_2 %>% as.data.frame() %>% 
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



bseq_varPart["ENSG00000116273",]
