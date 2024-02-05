# build heatmap and tree for Fo effectors

# install.packages("ggtree")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree")
# install.packages('devtools')
# devtools::install_github("YuLab-SMU/ggtreeExtra")
# install.packages('tidytree')

# Load required libraries
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(ape)
library(dplyr)
library(tidyr)
library(phytools)
library(tidytree)
library(textshape)
library(RColorBrewer)
library(ggnewscale)

# ---- Set dirs and load files ---- #
setwd("/Volumes/Jamie_EXT/Projects/Maei/exp/AllFusAnalysis-Chap3/R/")
# Specify  path to  Newick file
tree_file <- "./MaeiTEFPhylo.treefile"
# Load heatmap data matrix from the CSV file
data <- read.csv("/PhyloHeatmapData.csv")  # Adjust the path accordingly
# Load phylogeny meta date from the CSV file
metadata <- read.csv("./PhyloMetaData.csv") 

# ---- prepare metadata ---- #
# ensure its a df
metadata <- as.data.frame(metadata)
# clear empty values
metadata$fsp[metadata$fsp==""] <- NA # set empty cells in fsp column to NA
metadata$race[metadata$race==""] <- NA # set empty cells in race column to NA
# build df for just race 
race_df <- data.frame("race" = metadata[,c("race")])
rownames(race_df) <- metadata$label
# build df for fsp. 
fsp_df <- data.frame("fsp" = metadata[,c("fsp")])
rownames(fsp_df) <- metadata$label

# ---- Prepare the tree ---- #
# Read the phylogenetic tree from the Newick file
unrootedtree <- read.tree(tree_file)
# root the tree
tree <- root(unrootedtree, outgroup = c("F._graminearum_PH-1"))
# Build tree skeleton
p <- ggtree(tree,ladderize = T)  %<+% metadata

# ---- View the tree ---- #
# Useful for visualusing nodes etc 
ggtree(tree) + geom_tiplab(aes(label=node), align = T)
# save as file for viewing later. 
ggsave("outtree_nodes.png", width = 30, height = 15)

# ---- Prep the heatmap data ---- #
# remove the row names temporarily 
rownames_mat<- data[,1]
mat_data<- as.matrix(data[,-1])
#make data frame binary 
binary_matrix <- as.matrix(mat_data)
binary_matrix[binary_matrix > 0] <- 1
#put rownames back
rownames(binary_matrix)<-rownames_mat

# ---- Cluster the heatmap data ---- #
# Compute hierarchical clustering of columns
heatmap_dat <- cluster_matrix(binary_matrix, dim = 'col', method ="ward.D2")

# ---- Build the tree plot ---- #

p2 <- p + 
  geom_treescale(x = 0, y = 45, width = 0.004) + 
  geom_tiplab(aes(label = label), size = 6, offset = 0.005) +
  geom_tiplab(aes(label = race), color = "black", offset = 0.07, size = 6, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source), size = 5) +
  geom_rootedge()
p2

##save basic tree 
#ggsave("outtree_basic.png", width = 30, height = 15)

# ---- Build the full figure ---- #

# add effector heatmap
p3 <-gheatmap(p2, heatmap_dat, offset=0.125, colnames=FALSE, legend_title="Presence/Absence", color = "grey",  width = 2)  +
  scale_fill_continuous(name = "Presence/Absence",
                        low = "white", high = "lightblue",
                        breaks = c(0,1.00),
                        na.value = "grey")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

# add extra discrete scale
p4 <- p3 + new_scale_fill()
# add race data
p5 <- gheatmap(p4, fsp_df,
               offset = 0.062, 
               width = 0.03,
               color = "black",
               colnames = FALSE) +
  scale_fill_manual(name = "Fsp",
                    values = c("antiquewhite2","aquamarine2","azure2","coral2","chocolate2","darkgoldenrod2", "cornsilk", "lightblue", "darkolivegreen3", "darkslategrey", "lavender", "lightpink", "palegreen4", "tomato")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))
p5

# save final tree
ggsave("outtree.png", width = 30, height = 20)

# ---- Explore specific clades ---- #
#lactucae clade
p_lac <- viewClade(p2, node = 55)

gheatmap(p_lac, heatmap_dat, offset=0.06, colnames=FALSE, legend_title="Presence/Absence", color = NULL)  +
  theme(legend.position = "bottom")

#cubense clades


#dev.off()

