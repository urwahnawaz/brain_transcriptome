


## load data 

hca= read.csv("../../../BrainData/single-cell/hca/HCA-metadata_fixed.csv", header=TRUE)
colnames(hca)

vel = read.csv("../../../BrainData/single-cell/velmeshev/Velmeshev-metadata.csv", header= TRUE)
colnames(vel)

### hca data 

## summarise cell-types and also add percentages 

hca.CT = hca %>% group_by(subclass, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    mutate(percentage = (n/ sum(n)) * 100) 

hca.CT$MajorCelltype[hca.CT$MajorCelltype =="Astrocytes" | hca.CT$MajorCelltype == "Microglia" | 
                     hca.CT$MajorCelltype == "Oligodendrocytes" | hca.CT$MajorCelltype == "OPCs"] <- c("Neuroglia")

hca.CT %<>%
arrange(MajorCelltype, n) 


hca.CT$id <- seq(1, nrow(hca.CT))
label_data <- hca.CT

number_of_bar <- nrow(label_data)
number_of_bar
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
angle
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)




ggplot(hca.CT, aes(x=as.factor(id), y=n, fill=MajorCelltype)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", color = "white") +
    theme_minimal() + ylim(-100,12) +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) + scale_y_log10(limits = c(1e-2, NA)) +
    coord_polar() + 
    geom_text(data=label_data, 
              aes(x=id, y=n+10, label=subclass, hjust=hjust), 
              color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

### circular plot for velmeshev data 


vel.CT =vel %>%  group_by(CellType, MajorCelltype) %>% 
    summarise(n = n()) %>% 
    as.data.frame() %>% 
    arrange(MajorCelltype, n)
vel.CT
