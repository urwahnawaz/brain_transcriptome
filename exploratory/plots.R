## Figures 
libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl", 
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl", "corrplot", "viridis", 
         "sunburstR")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

col_palette = c("#162424","#476460", "#7f9690", "#b8c0a9","#88947e","#373e36", "#9dac81","#7e979e", "#c9cdce", "#a4b4b1", 
                "#8b8d8c","#ce9c83","#6d2010", "#9b4129","#dfb29b","#e0d7c8", "#f1e3e3","#e6c2a8","#d4ae87", "#f5c39e", 
                "#543320", "#532414","#b87356","#cb8b70","#dca996","#f4d6be","#d1a88a","#be825d","#d6a392","#966f4e")

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


### SN plots 
hca= read.csv("../../../BrainData/single-cell/hca/HCA-metadata_fixed.csv", header=TRUE)
vel = read.csv("../../../BrainData/single-cell/velmeshev/Velmeshev-metadata.csv", header= TRUE)

### circle bar graph 


# Create dataset
data <- data.frame(
    individual=paste( "Mister ", seq(1,60), sep=""),
    group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
    value=sample( seq(10,100), 60, replace=T)
)

data %<>% arrange(group, value)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))


# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


data

# Make the plot
# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", alpha=0.5) +
    ylim(-100,120) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p


### HCA 

hca.CT = hca %>% group_by(subclass, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    arrange(MajorCelltype, n)

hca.CT$FinalPath[hca.CT$MajorCelltype =="Astrocytes" | hca.CT$MajorCelltype == "Microglia" | 
                         hca.CT$MajorCelltype == "Oligodendrocytes" | hca.CT$MajorCelltype == "OPCs"] <- c("Neuroglia")


hca.CT %>% 
    mutate(path = paste(FinalPath,MajorCelltype, subclass, sep = "-")) %>% 
    dplyr::select(path, n) %>%
  #  mutate_at(.vars ="path", .funs = gsub, pattern = "\\^- -", replacement = "-")
    sunburst(., legend = FALSE)

hca.CT$FinalPath[grepl("Neurons", hca.CT$MajorCelltype, ignore.case = TRUE)] <- c("Neurons")
hca.CT$FinalPath[grepl("Vasculature", hca.CT$MajorCelltype, ignore.case = TRUE)] <- c("Vasculature")
hca.CT$MajorCelltype[hca.CT$MajorCelltype =="Astrocytes" | hca.CT$MajorCelltype == "Microglia" | 
                     hca.CT$MajorCelltype == "Oligodendrocytes" |
                         hca.CT$MajorCelltype == "OPCs" | hca.CT$MajorCelltype == "Vasculature"] <- NA



hca.CT
hca.CT$FinalPath[is.na(hca.CT$FinalPath)] <- " "



table(vel.CT$MajorCelltype) %>% length()

empty_bar <- 7
to_add <- data.frame( matrix(NA, empty_bar*nlevels(vel.CT$group), ncol(vel.CT)) )
colnames(to_add) <- colnames(vel.CT)
to_add$MajorCelltype <- rep(levels(vel.CT$MajorCelltype), each=empty_bar)
vel.CT <- rbind(vel.CT, to_add)
vel.CT <- vel.CT %>% arrange(MajorCelltype)
vel.CT$id <- seq(1, nrow(vel.CT))


label_data <- vel.CT
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


ggplot(vel.CT, aes(x=as.factor(id), y=n, fill=MajorCelltype)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", alpha=0.5) +
    ylim(-100,100) +
    theme_minimal() + scale_y_log10() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
   geom_text(data=label_data, 
              aes(x=id, y=n+10, label=subclass, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 



vel.CT =vel %>%  group_by(CellType, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    arrange(MajorCelltype, n)




empty_bar <- 9
to_add <- data.frame( matrix(NA, empty_bar*nlevels(vel.CT$group), ncol(vel.CT)) )
colnames(to_add) <- colnames(vel.CT)
to_add$MajorCelltype <- rep(levels(vel.CT$MajorCelltype), each=empty_bar)
vel.CT <- rbind(vel.CT, to_add)
vel.CT <- vel.CT %>% arrange(MajorCelltype)
vel.CT$id <- seq(1, nrow(vel.CT))


label_data <- vel.CT
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

## convert the values to percentages?


ggplot(vel.CT, aes(x=as.factor(id), y=n, fill=MajorCelltype)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", alpha=0.5) +
    ylim(-10,6000) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,10), "cm") 
    ) +
    coord_polar() + 
    geom_text(data=label_data, 
              aes(x=id, y=n+10, label=CellType, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 


## sunburstR plots

data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/11_SevCatOneNumNestedOneObsPerGroup.csv", header=T, sep=";")
head(data)
data[ which(data$value==-1),"value"] <- 1
colnames(data) <- c("Continent", "Region", "Country", "Pop")

# Reformat data for the sunburstR package
data <- data %>%
    filter(Continent != "") %>%
    mutate(path = paste(Continent, Region, Country, sep="-")) %>%
    dplyr::select(path, Pop)
data



vel.CT %>% 
    mutate_at(.vars = "CellType", .funs=gsub, pattern = "\\-", replacement = " ") %>% 
    mutate(path = paste( MajorCelltype, CellType, sep = "-")) %>% 
    dplyr::select(path, n) %>% 
    sunburst(., legend = FALSE)

# Plot

## make one subburst plot 
p <- sunburst(data, legend=FALSE)
p
?sunburst

sum(hca.CT$n)
?colSums
hca.CT %<>% 
    mutate(percentage = (n/sum(n) * 100)) 


sum(hca.CT$percentage)

library(dplyr)
eg= iris %>% 
    slice(1:4) #%>% 
    mutate(Test=Sepal.Length/45,Test=scales::percent(Test)) 

iris

