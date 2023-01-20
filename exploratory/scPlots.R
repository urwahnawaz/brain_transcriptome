# single-cell RNA-seq plots 
HCA.md = read.csv("../../Data/SN Formatted /HCA-metadata_fixed.csv", header= TRUE)
head(HCA.md)
table(HCA.md$cluster)
HCA.md %>% 
  group_by(MajorCelltype) %>%
  summarise(counts = n()) %>% 
  ggplot(aes(x = MajorCelltype, y = counts)) +
  geom_bar(fill = "#88738D", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) + coord_flip() +
  theme_pubclean() + ylab("# Cells") + xlab("Major Cell-Types")
 table(HCA.md$AgeNumeric)

table(HCA.md$full_genotype)
HCA.md
Vel.md = read.csv("../../Data/SN Formatted /Velmeshev-metadata.csv", header= TRUE)
table(Vel.md$DonorID, Vel.md$AgeNumeric)
unique(Vel.md$DonorID) ## 31 donors
Vel.md %>%
  group_by(MajorCelltype) %>%
  summarise(counts = n()) %>% 
  ggplot(aes(x = MajorCelltype, y = counts)) +
  geom_bar(fill = "#88738D", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) + coord_flip() +
  theme_pubclean() + ylab("# Cells") + xlab("Major Cell-Types")
