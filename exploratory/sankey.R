## Sankey plots for total regions and structures in the dataset 
library(networkD3)
library(ggsankey)
library(ggplot2)
library(dplyr)
ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)
## what a sankey df looks like 
nodes = data.frame("name" = 
                     c("Node A", # Node 0
                       "Node B", # Node 1
                       "Node C", # Node 2
                       "Node D"))# Node 3
nodes

links = as.data.frame(matrix(c(
  0, 1, 10, # Each row represents a link. The first number
  0, 2, 20, # represents the node being conntected from. 
  1, 3, 30, # the second number represents the node connected to.
  2, 3, 40),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
names(links) = c("source", "target", "value")
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30)
links


##

### 
## gtex data 
gtex.md = read.csv("../../Data/FormattedData/GTEx-metadata.csv", header=TRUE) %>% 
  dplyr::select(StructureAcronym, Regions) %>% 
  mutate(dataset = "GTEx")
bspan.md = read.csv("../../Data/FormattedData/BrainSpan-metadata.csv", header= TRUE) %>% 
  dplyr::select(StructureAcronym, Regions) %>% 
  mutate(dataset = "BrainSpan")

bseq.md= read.csv("../../Data/FormattedData/BrainSeq-metadata.csv", header=TRUE) %>% 
  dplyr::select(StructureAcronym, Regions) %>% 
  mutate(dataset = "BrainSeq")
regions_sankey = rbind(gtex.md, bspan.md, bseq.md)

table(gtex.md$StructureAcronym, gtex.md$Regions)

regions_sankey 


regions_sankey %>% make_long( StructureAcronym, Regions) %>% 
  ggplot(aes(x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = factor(node), 
             label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1) +
  theme_sankey(base_size = 16) +
  theme(legend.position= "none")
