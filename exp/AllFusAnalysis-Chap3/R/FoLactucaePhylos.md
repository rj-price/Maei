Fola Six Gene Trees
================
JamiePike
2024-02-16

## Build SIX gene trees for Fola SIX gene phylogenies

Phylogenies were generated as part of the FoLac project with NIAB.

``` r
# ---- Set dirs and load files ---- #
setwd("/Volumes/Jamie_EXT/Projects/Maei/exp/AllFusAnalysis-Chap3/R/")
# Specify  path to  Newick file
SIX8_tree_file <- "./Fola_SIX8_Phylo1.treefile"
SIX9_tree_file <-  "./Fola_SIX9_Phylo1.treefile"
SIX14_tree_file <- "./Fola_SIX14_Phylo1.treefile"
# Load phylogeny meta date from the CSV file
metadata <- read.csv("./PhyloMetaData.csv") 
```

After loading my data, I prepared it for subsequent analysis. The
metadata and heatmap data are currently in .csv format, but that saves
empty cells as ““, and I want them to be `NA` values (particularly the
metadata). I also added columns which contains the”full name” for each
genome assembly, where full name is the genus (abbreviated), species,
and forma speciales (fsp).

``` r
# ---- prepare metadata ---- #
# ensure its a df
metadata <- as.data.frame(metadata)
# clear empty values
metadata$fsp[metadata$fsp==""] <- NA # set empty cells in fsp column to NA
metadata$race[metadata$race==""] <- NA # set empty cells in race column to NA

#create a column for full fsp.
metadata <- metadata %>% unite("full_name", c(species,fsp), sep = " fsp. ", remove = F, na.rm = T)

#create a column for full fsp and isolate code
metadata <- metadata %>% unite("full_ID", c(full_name,isolate_code), sep = " ", remove = F, na.rm = T)

# ---- add reference sequences ---- #

# add in the reference sequences
metadata <- metadata %>% 
  add_row(label = "Fo._fsp._lycopersici_4287_SIX8_Reference_FJ755837.1", full_ID = "F. oxysporum fsp. lycopersici 4287 SIX8 Reference FJ755837.1", full_name = "F. oxysporum fsp. lycopersici", fsp = "lycopersici", race = "Race 2", .before = 1) 

metadata <- metadata %>% 
  add_row(label = "Fo._fsp._lycopersici_4287_SIX9_Reference_KC701447.1", full_ID = "F. oxysporum fsp. lycopersici 4287 SIX9 Reference KC701447.1", full_name = "F. oxysporum fsp. lycopersici", fsp = "lycopersici", race = "Race 2", .before = 2)

metadata <- metadata %>% 
  add_row(label = "Fo._fsp._lycopersici_4287_SIX14_Reference_KC701452.1", full_ID = "F. oxysporum fsp. lycopersici 4287 SIX14 Reference KC701452.1", full_name = "F. oxysporum fsp. lycopersici", fsp = "lycopersici", race = "Race 2", .before = 2)


# build df for fsp so we can add it as a colour scale to the tree
fsp_df <- data.frame("full_name" = metadata[,c("full_name")] )
rownames(fsp_df) <- metadata$label
```

## SIX8 tree

``` r
# ---- Prepare the tree ---- #
# Read the phylogenetic tree from the Newick file
SIX8_unrootedtree <- read.newick(SIX8_tree_file, node.label='label' )
# root the tree
SIX8_tree <- root(SIX8_unrootedtree, outgroup = c("Fo._fsp._lycopersici_4287_SIX8_Reference_FJ755837.1"))

#lets drop the tips from the fsp not indcluded in the maie analysis but included in the NIAB project
# first make the tree data a data frame
SIX8_tree_df <- SIX8_tree %>% as_tibble()

# keep only the species names for now
SIX8_tree_df1 <- SIX8_tree_df %>% 
  filter(grepl("^Fo\\.[^0-9]", label))

# now find the differences between the two data frames
tips_to_remove <- setdiff(SIX8_tree_df1$label, metadata$label)
# prune teh tree
SIX8_tree_pruned <- drop.tip(SIX8_tree, tips_to_remove)

# ---- Build tree skeleton ---- #
p <- ggtree(SIX8_tree_pruned, ladderize = F)  %<+% metadata


# ---- View the tree ---- #
# Useful for visualusing nodes etc 
p_nodes <- p + 
  geom_point(data = td_filter(isTip %in% metadata$label)) +
   geom_tiplab(aes(label = label), offset = 0.005) + 
  geom_text2(aes(label = parent), hjust = -0.1, size = 3)+ # add node names
  coord_cartesian(clip = "off") # stop names being trimmed off
```

Now I have my basic tree, I can start to build something that will stand
alone. First I added the metadata (full name, the isolate code, and
race).

``` r
# ---- Build the tree plot ---- #

p2 <- p + 
  geom_treescale(x = 0, y = -0.5, width = 0.1) + 
  geom_tiplab(aes(label = full_ID), offset = 0.005) +
  geom_tiplab(aes(label = race), color = "grey20", offset = 0.2, linetype = "blank", geom = "text", align = TRUE, hjust = 1)+
  geom_tippoint(aes(shape = source), size = 2.5) +
  geom_rootedge() +
  theme(legend.position = "bottom" ) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 60), nudge_x = -0.0012, nudge_y = -0.25) 

# Add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()

# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.2, 
               width = 0.03,
               color = "grey20",
               colnames = FALSE) +
  scale_fill_manual(name = "Forma specialis",
                    values = c( "goldenrod","gold", "tan","darkolivegreen3",  "tomato", "lavender"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 5))) 
```

Now plot and save it.

``` r
# plot it 
plot(p4)
```

<figure>
<img src="FoLactucaePhylos_files/figure-gfm/TEF1a%20tree%20plot-1.png"
alt="SIX8 phylogeny of Fusarium isolates included int effector anaylisis.The tree is rooted through Foly reference." />
<figcaption aria-hidden="true">SIX8 phylogeny of Fusarium isolates
included int effector anaylisis.The tree is rooted through Foly
reference.</figcaption>
</figure>

``` r
##save basic tree 
ggsave("lactucaeSIX8tree.png", width = 10, height = 15)
```

## SIX9 tree

``` r
# ---- Prepare the tree ---- #
# Read the phylogenetic tree from the Newick file
SIX9_unrootedtree <- read.newick(SIX9_tree_file, node.label='label' )
# root the tree
SIX9_tree <- root(SIX9_unrootedtree, outgroup = c("Fo._fsp._lycopersici_4287_SIX9_Reference_KC701447.1"))

#lets drop the tips from the fsp not indcluded in the maie analysis but included in the NIAB project
# first make the tree data a data frame
SIX9_tree_df <- SIX9_tree %>% as_tibble()

# keep only the species names for now
SIX9_tree_df1 <- SIX9_tree_df %>% 
  filter(grepl("^Fo\\.[^0-9]", label))

# now find the differences between the two data frames
tips_to_remove <- setdiff(SIX9_tree_df1$label, metadata$label)
# prune teh tree
SIX9_tree_pruned <- drop.tip(SIX9_tree, tips_to_remove)

# ---- Build tree skeleton ---- #
p <- ggtree(SIX9_tree_pruned, ladderize = F)  %<+% metadata


# ---- View the tree ---- #
# Useful for visualusing nodes etc 
p_nodes <- p + 
  geom_tiplab(aes(label = full_name), offset = 0.05) + 
  geom_text2(aes(label = node), hjust = -0.1, size = 3, offeset =0.02)+ # add node names
  coord_cartesian(clip = "off") # stop names being trimmed off
```

Now I have my basic tree, I can start to build something that will stand
alone. First I added the metadata (full name, the isolate code, and
race).

``` r
# ---- Build the tree plot ---- #

p2 <- p + 
  geom_cladelabel(node = 52, label = "Group 1", offset = 0.55) +
  geom_cladelabel(node = 64, label = "Group 2", offset = 0.55) +
  geom_cladelabel(node = 76, label = "Group 3", offset = 0.55) +
  geom_cladelabel(node = 90, label = "Group 4", offset = 0.55) +
  geom_treescale(x = 0, y = -0.5, width = 0.1) + 
  geom_tiplab(aes(label = full_ID), offset = 0.005) +
  geom_tiplab(aes(label = race), color = "grey20", offset = 0.75, linetype = "blank", geom = "text", align = TRUE, hjust = 1)+
  geom_tippoint(aes(shape = source), size = 2.5) +
  geom_rootedge() +
  theme(legend.position = "bottom" ) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 60), nudge_x = -0.05) 

# Add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()

# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.75, 
               width = 0.03,
               color = "grey20",
               colnames = FALSE) +
  scale_fill_manual(name = "Forma specialis",
                    values = c("blue","goldenrod","grey90","gold", "tan","darkolivegreen3", "indianred", "tomato", "lavender", "yellow", "palegreen4", "slateblue", "steelblue"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 5))) 
```

Now plot and save it.

``` r
# plot it 
plot(p4)
```

<figure>
<img src="FoLactucaePhylos_files/figure-gfm/SIX9%20tree%20plot-1.png"
alt="SIX9 phylogeny of Fusarium isolates included int effector anaylisis.The tree is rooted through Foly reference." />
<figcaption aria-hidden="true">SIX9 phylogeny of Fusarium isolates
included int effector anaylisis.The tree is rooted through Foly
reference.</figcaption>
</figure>

``` r
##save basic tree 
ggsave("lactucaeSIX9tree.png", width = 12, height = 20)
```

## SIX14 tree

``` r
# ---- Prepare the tree ---- #
# Read the phylogenetic tree from the Newick file
SIX14_unrootedtree <- read.newick(SIX14_tree_file, node.label='label' )
# root the tree
SIX14_tree <- root(SIX14_unrootedtree, outgroup = c("Fo._fsp._lycopersici_4287_SIX14_Reference_KC701452.1"))

#lets drop the tips from the fsp not indcluded in the maie analysis but included in the NIAB project
# first make the tree data a data frame
SIX14_tree_df <- SIX14_tree %>% as_tibble()

# keep only the species names for now
SIX14_tree_df1 <- SIX14_tree_df %>% 
  filter(grepl("^Fo\\.[^0-9]", label))

# now find the differences between the two data frames
tips_to_remove <- setdiff(SIX14_tree_df1$label, metadata$label)
# prune teh tree
SIX14_tree_pruned <- drop.tip(SIX14_tree, tips_to_remove)

# ---- Build tree skeleton ---- #
p <- ggtree(SIX14_tree_pruned, ladderize = F)  %<+% metadata


# ---- View the tree ---- #
# Useful for visualusing nodes etc 
p_nodes <- p + 
  geom_tiplab(aes(label = full_name), offset = 0.05) + 
  geom_text2(aes(label = node), hjust = -0.1, size = 3, offeset =0.02)+ # add node names
  coord_cartesian(clip = "off") # stop names being trimmed off
```

``` r
# ---- Build the tree plot ---- #

p2 <- p + 
  geom_cladelabel(node = 29, label = "Group 1", offset = 0.55) +
  geom_cladelabel(node = 37, label = "Group 2", offset = 0.55) +
  geom_treescale(x = 0, y = -0.5, width = 0.1) + 
  geom_tiplab(aes(label = full_ID), offset = 0.005) +
  geom_tiplab(aes(label = race), color = "grey20", offset = 0.75, linetype = "blank", geom = "text", align = TRUE, hjust = 1)+
  geom_tippoint(aes(shape = source), size = 2.5) +
  geom_rootedge() +
  theme(legend.position = "bottom" ) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 60), nudge_x = -0.05) 

# Add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()

# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.75, 
               width = 0.03,
               color = "grey20",
               colnames = FALSE) +
  scale_fill_manual(name = "Forma specialis",
                    values = c("blue","grey90", "tan","darkolivegreen3", "indianred", "tomato", "slateblue"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 5))) 
```

Now plot it

``` r
# plot it 
plot(p4)
```

<figure>
<img src="FoLactucaePhylos_files/figure-gfm/SIX14%20tree%20plot-1.png"
alt="SIX14 phylogeny of Fusarium isolates included int effector anaylisis.The tree is rooted through Foly reference." />
<figcaption aria-hidden="true">SIX14 phylogeny of Fusarium isolates
included int effector anaylisis.The tree is rooted through Foly
reference.</figcaption>
</figure>

``` r
##save basic tree 
ggsave("lactucaeSIX14tree.png", width = 12, height = 20)
```
