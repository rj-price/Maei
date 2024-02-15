Analysis of Candidate Effectors
================

## Data set

Candidate effectors were generated as part of the Third Chapter of my
thesis. They were identified using the [Maei
pipeline](https://github.com/JamiePike/Maei) in a set of *Fusarium*
genomes. For the full set of genomes and the previous command line
scripts and analysis, see
[Process.md](https://github.com/JamiePike/Maei/docs/Process.md).

``` r
# ---- Set dirs and load files ---- #
setwd("/Volumes/Jamie_EXT/Projects/Maei/exp/AllFusAnalysis-Chap3/R/")
# Specify  path to  Newick file
tree_file <- "MaeiTEFPhylo.phylo2.treefile"
# Load heatmap data matrix from the CSV file
data <- read.csv("./PhyloHeatmapData.csv")  # Adjust the path accordingly
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


# build df for fsp so we can add it as a colour scale to the tree
fsp_df <- data.frame("full_name" = metadata[,c("full_name")] )
rownames(fsp_df) <- metadata$label
```

## Tef-1a tree of isolates

Now I have my database of isolates, I need to know how they are related.
I have already performed the phylogenetic analysis, here I am plotting
the results. First, I rooted the tree through F.\_graminearum_PH-1
*Tef1-a*, and built the basic tree skeleton. I also loaded in the
metadata for labeling later. Additionally, I built a basic tree with the
nodes labeled, so that I can use that for later reference

``` r
# ---- Prepare the tree ---- #
# Read the phylogenetic tree from the Newick file
unrootedtree <- read.newick(tree_file, node.label='label' )
# root the tree
tree <- root(unrootedtree, outgroup = c("F._graminearum_PH-1"))

# Identify the node number of "F._graminearum_PH-1"
node_number <- which(tree$tip.label == "F._graminearum_PH-1")
# Trim the F._graminearum_PH-1 branch
m <- MRCA(tree, 48)
tree_with_group <- groupClade(tree, m)

# ---- Build tree skeleton ---- #
p <- ggtree(tree, ladderize = T)  %<+% metadata

# Adjust the branch length for "F._graminearum_PH-1"
p$data[p$data$node %in% node_number, "x"] <- mean(p$data$x)

# ---- View the tree ---- #
# Useful for visualusing nodes etc 
p_nodes <- p + 
  geom_text2(aes(label = parent), hjust = -0.1, size = 3)+ # add node names
  geom_tiplab(aes(label = label), offset = 0.005) +
  geom_cladelabel(node=node_number, label="Adjusted by: 0.03", 
                  color='black') +
  coord_cartesian(clip = "off") # stop names being trimmed off
```

Now I have my basic tree, I can start to build something that will stand
alone. First I added the metadata (full name, the isolate code, and
race).

``` r
# ---- Build the tree plot ---- #

p1 <- p +
  xlim(0,0.045)

p2 <- p + 
  geom_highlight(node = 88, fill = "mistyrose") + # colour the Fs node
  geom_highlight(node= 49, fill = "lemonchiffon1" ) + # colour the fo node
  geom_treescale(x = 0, y = 1, width = 0.004) + 
  geom_tiplab(aes(label = full_ID), offset = 0.0004) +
 # geom_tiplab(aes(label = isolate_code), color = "grey20", offset = 0.0045, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), color = "grey20", offset = 0.012, linetype = "blank", geom = "text", align = TRUE, hjust = 1)+
  geom_tippoint(aes(shape = source), size = 2.5) +
  geom_rootedge() +
  theme(legend.position = "bottom" ) +
  geom_cladelabel(node=node_number, label="Adjusted by: 0.03", offset = 0.008, color='black') + geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 60), nudge_x = -0.0012, nudge_y = -0.25) 

# Add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()

# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.012, 
               width = 0.03,
               color = "grey20",
               colnames = FALSE) +
  scale_fill_manual(name = "Forma specialis",
                    values = c( "burlywood4", "grey","blue","purple","goldenrod","grey90","gold","brown", "tan","darkolivegreen3", "indianred", "tomato", "lavender", "yellow", "palegreen4", "slateblue", "steelblue", "pink", "red"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 5))) 
```

The TR4 isolates group into a monophyletic clade, as expected; as do the
*lactucae* isolates. The TNAU isolates group as they did in the previous
TEF1a and RBP2 phylogenies (chapter 2). Isolates of the same fsp and
race do not all group in the same clade (e.g. *Fo.* fsp. *cubense* Race
1)

``` r
# plot it 
plot(p4)
```

<figure>
<img
src="AnalysisCandidateEffectorSets_files/figure-gfm/TEF1a%20tree%20plot-1.png"
alt="TEF phylogeny of Fusarium isolates included int effector anaylisis.The tree is rooted through F. graminearum." />
<figcaption aria-hidden="true">TEF phylogeny of Fusarium isolates
included int effector anaylisis.The tree is rooted through F.
graminearum.</figcaption>
</figure>

``` r
##save basic tree 
ggsave("BasicTEFPhylo.png", width = 12, height = 20)
```

## Analysis of candidate effectors, *mimps* and genome size in data set

Once the candidate effectors had been identified in the candidates, I
conducted some basic statistics to see how they were distributed.

First, I subset the data I needed from the `Metadata`, and removed the
isolates that were included in the TEF1a phylogenies but not the Maei
analysis.

``` r
# ----  tidy the data for stats ----#
#subset the df so that we can perform stats on genome size and mimp/cand eff distribution
stats_data <- select(metadata, "species", "fsp", "isolate_code","no._mimps","no._cand_effs","genome_size") %>%
#Rename the columns to reduce the long titles
   dplyr::rename(isolate=isolate_code,
         mimps=no._mimps,
         candidate_effectors=no._cand_effs,
         assembly_size =genome_size)

#we need to drop rows which were not included in the Maei analysis
stats_data <- stats_data %>%
  drop_na(candidate_effectors)

# view the table for stats data
knitr::kable(stats_data)
```

| species        | fsp          | isolate      | mimps | candidate_effectors | assembly_size |
|:---------------|:-------------|:-------------|------:|--------------------:|--------------:|
| F. oxysporum   | apii         | 207A         |   420 |                 388 |          64.6 |
| F. oxysporum   | apii         | 274AC        |   217 |                 357 |          67.3 |
| F. oxysporum   | apii         | AJ498        |   539 |                 332 |          64.6 |
| F. oxysporum   | apii         | AJ720        |   442 |                 399 |          64.7 |
| F. oxysporum   | apii         | NRRL38295    |   210 |                 328 |          65.3 |
| F. oxysporum   | cepae        | FoC_Fus2     |   325 |                 359 |          53.4 |
| F. oxysporum   | conglutinans | Fo5176       |   442 |                 385 |          68.0 |
| F. oxysporum   | coriandrii   | 3-2          |   478 |                 315 |          65.4 |
| F. oxysporum   | coriandrii   | AJ615        |   675 |                 603 |          69.3 |
| F. oxysporum   | cubense      | 58           |   165 |                 108 |          48.2 |
| F. oxysporum   | cubense      | 60           |    39 |                  45 |          48.6 |
| F. oxysporum   | cubense      | 160527       |   183 |                  95 |          51.1 |
| F. oxysporum   | cubense      | B2           |   142 |                  60 |          52.9 |
| F. oxysporum   | cubense      | C1HIR_9889   |   145 |                  70 |          46.7 |
| F. oxysporum   | cubense      | N2           |   141 |                  63 |          47.7 |
| F. oxysporum   | cubense      | NRRL_54006   |   105 |                  70 |          46.6 |
| F. oxysporum   | cubense      | Pers4        |    87 |                  62 |          46.4 |
| F. oxysporum   | cubense      | UK0001       |   160 |                 127 |          48.6 |
| F. oxysporum   | cubense      | VPRI44079    |   149 |                 107 |          49.5 |
| F. oxysporum   | cubense      | VPRI44081    |   137 |                  45 |          47.2 |
| F. oxysporum   | cubense      | VPRI44082    |   146 |                  50 |          46.3 |
| F. oxysporum   | cubense      | VPRI44083    |   146 |                  50 |          46.3 |
| F. oxysporum   | cubense      | VPRI44084    |   179 |                 104 |          49.5 |
| F. oxysporum   | endophyte    | Fo47         |    70 |                  81 |          50.4 |
| F. oxysporum   | from rocket  | AJ174        |   419 |                 169 |          62.6 |
| F. oxysporum   | lactucae     | AJ516        |   522 |                 548 |          68.8 |
| F. oxysporum   | lactucae     | AJ520        |   533 |                 343 |          62.2 |
| F. oxysporum   | lactucae     | AJ592        |   615 |                 482 |          66.0 |
| F. oxysporum   | lactucae     | AJ705        |   614 |                 490 |          66.2 |
| F. oxysporum   | lactucae     | AJ718        |   536 |                 260 |          62.1 |
| F. oxysporum   | lactucae     | AJ865        |   569 |                 296 |          62.7 |
| F. oxysporum   | lini         | 39_C0058     |   263 |                 237 |          59.2 |
| F. oxysporum   | lycopersici  | 4287         |   332 |                 337 |          56.2 |
| F. oxysporum   | matthiolae   | AJ260        |   301 |                 209 |          60.3 |
| F. oxysporum   | narcissus    | FON63        |   555 |                 251 |          60.0 |
| F. oxysporum   | niveum       | 110407-3-1-1 |    52 |                  39 |          49.7 |
| F. oxysporum   | rapae        | Tf1208       |   377 |                 278 |          59.8 |
| F. oxysporum   | vasinfectum  | TF1          |   199 |                  91 |          50.0 |
| F. sacchari    | NA           | FS66         |     9 |                  12 |          47.5 |
| F. sacchari    | NA           | NRRL_66326   |    10 |                  15 |          42.8 |
| F. graminearum | NA           | PH-1         |     1 |                   6 |          38.0 |
| Fusarium       | NA           | SY-2         |     3 |                  12 |          44.2 |

The number of *mimps*, candidate effectors identified, and genome sizes
varies between genome assemblies.

``` r
# ---- sumarise the data ---- #
summary(stats_data$mimps)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     1.0   141.2   204.5   277.4   442.0   675.0

``` r
summary(stats_data$candidate_effectors)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    6.00   62.25  148.00  206.62  335.75  603.00

``` r
summary(stats_data$assembly_size)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   38.00   47.83   53.15   55.40   64.12   69.30

I also want to look at specific fsp. to determine if the avernage number
of effectors, *mimps*, and genome size is smaller compared to other fsp.
To do this, I `subset` the fsp. I am interested in and then use the
`summary` function in R.

``` r
# ---- subset the fsp of interest ---- #
# cubense
foc_stats <- subset(metadata, grepl("Fo._fsp._cubense", label))
summary(foc_stats)
```

    ##     label             full_ID           full_name           species         
    ##  Length:14          Length:14          Length:14          Length:14         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  species_group          fsp            isolate_code           race          
    ##  Length:14          Length:14          Length:14          Length:14         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     source            no._mimps     no._cand_effs     genome_size   
    ##  Length:14          Min.   : 39.0   Min.   : 45.00   Min.   :46.30  
    ##  Class :character   1st Qu.:138.0   1st Qu.: 52.50   1st Qu.:46.62  
    ##  Mode  :character   Median :145.5   Median : 66.50   Median :47.95  
    ##                     Mean   :137.4   Mean   : 75.43   Mean   :48.26  
    ##                     3rd Qu.:157.2   3rd Qu.:101.75   3rd Qu.:49.27  
    ##                     Max.   :183.0   Max.   :127.00   Max.   :52.90

``` r
# lactucae
fola_stats <- subset(metadata, grepl("Fo._fsp._lactucae", label))
summary(fola_stats)
```

    ##     label             full_ID           full_name           species         
    ##  Length:6           Length:6           Length:6           Length:6          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  species_group          fsp            isolate_code           race          
    ##  Length:6           Length:6           Length:6           Length:6          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     source            no._mimps     no._cand_effs    genome_size   
    ##  Length:6           Min.   :522.0   Min.   :260.0   Min.   :62.10  
    ##  Class :character   1st Qu.:533.8   1st Qu.:307.8   1st Qu.:62.33  
    ##  Mode  :character   Median :552.5   Median :412.5   Median :64.35  
    ##                     Mean   :564.8   Mean   :403.2   Mean   :64.67  
    ##                     3rd Qu.:602.8   3rd Qu.:488.0   3rd Qu.:66.15  
    ##                     Max.   :615.0   Max.   :548.0   Max.   :68.80

``` r
# apii 
foa_stats <- subset(metadata, grepl("Fo._fsp._apii", label))
summary(foa_stats)
```

    ##     label             full_ID           full_name           species         
    ##  Length:5           Length:5           Length:5           Length:5          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  species_group          fsp            isolate_code           race          
    ##  Length:5           Length:5           Length:5           Length:5          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     source            no._mimps     no._cand_effs    genome_size  
    ##  Length:5           Min.   :210.0   Min.   :328.0   Min.   :64.6  
    ##  Class :character   1st Qu.:217.0   1st Qu.:332.0   1st Qu.:64.6  
    ##  Mode  :character   Median :420.0   Median :357.0   Median :64.7  
    ##                     Mean   :365.6   Mean   :360.8   Mean   :65.3  
    ##                     3rd Qu.:442.0   3rd Qu.:388.0   3rd Qu.:65.3  
    ##                     Max.   :539.0   Max.   :399.0   Max.   :67.3

``` r
# coriandrii
foci_stats <- subset(metadata, grepl("Fo._fsp._coriandrii", label)) 
foci_stats <- foci_stats %>%
  drop_na(no._mimps)
summary(foci_stats)
```

    ##     label             full_ID           full_name           species         
    ##  Length:2           Length:2           Length:2           Length:2          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  species_group          fsp            isolate_code           race          
    ##  Length:2           Length:2           Length:2           Length:2          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     source            no._mimps     no._cand_effs  genome_size   
    ##  Length:2           Min.   :478.0   Min.   :315   Min.   :65.40  
    ##  Class :character   1st Qu.:527.2   1st Qu.:387   1st Qu.:66.38  
    ##  Mode  :character   Median :576.5   Median :459   Median :67.35  
    ##                     Mean   :576.5   Mean   :459   Mean   :67.35  
    ##                     3rd Qu.:625.8   3rd Qu.:531   3rd Qu.:68.33  
    ##                     Max.   :675.0   Max.   :603   Max.   :69.30

### Correlation and distribution

I want to see if there is a relationship between the assembly size,
*mimp* content, and candidate effector number, so calculated
correlation.

First, I checked the distribution of the data, to determine the best
correlation test to perform.

``` r
# ---- Build the histograms ---- #
# visualise the mimp distribution
mimps_histo <- ggplot(stats_data, aes(x = mimps)) +
  geom_histogram(fill = "#0c4c8a", colour = "grey20" ) +
  theme_bw()

# visualise the mimp distribution
cands_histo <- ggplot(stats_data, aes(x = candidate_effectors)) +
  geom_histogram(fill = "#0c4c8a", colour = "grey20" ) +
  theme_bw()

# visualise the mimp distribution
size_histo <- ggplot(stats_data, aes(x = assembly_size)) +
  geom_histogram(fill = "#0c4c8a", colour = "grey20" ) +
  theme_bw()
```

The data are not normality distributed. The Pearson correlation is
computed by default with the `cor.test()` function.

I also visualized the data using scatter plots to check that the data
are linearly related.

``` r
# ---- Build the scatter plots ---- #
#Correlation between the total number of effectors and  the total number of mimps
#Visualize the relationship
effectors_v_mimps_relat <- ggplot(stats_data) +
  aes(x = candidate_effectors, y = mimps) +
  geom_point(colour = "#0c4c8a") +
  geom_smooth(method = "lm", se = FALSE, color = "grey20") +
  theme_bw()

#Correlation between  the total number of effectors and Assembly Size?
#Visualize the relationship
effectors_v_assembly_size_relat <- ggplot(stats_data) +
  aes(x = candidate_effectors, y = assembly_size) +
  geom_point(colour = "#0c4c8a") +
  geom_smooth(method = "lm", se = FALSE, color = "grey20") +
  theme_bw()

#Correlation between the total number of mimps and Assembly Size?
#Visualize the relationship
mimps_v_assembly_size_relat <- ggplot(stats_data) +
  aes(x = mimps, y = assembly_size) +
  geom_point(colour = "#0c4c8a") +
  geom_smooth(method = "lm", se = FALSE, color = "grey20") +
  theme_bw()
```

``` r
# combine the plots 
sats_plots <- ggarrange(mimps_histo, cands_histo, size_histo, mimps_v_assembly_size_relat, effectors_v_assembly_size_relat, effectors_v_mimps_relat,
          ncol = 3, 
          nrow = 2)

# plot the combined plots
plot(sats_plots)
```

<figure>
<img
src="AnalysisCandidateEffectorSets_files/figure-gfm/stats%20plots-1.png"
alt="Plots of meta data, showing distribution and relationship of assembly size, mimp, and candidate effector content" />
<figcaption aria-hidden="true">Plots of meta data, showing distribution
and relationship of assembly size, mimp, and candidate effector
content</figcaption>
</figure>

``` r
#save it 
ggsave("./SummaryStats.png", width = 20, height = 15)
```

There is a general positive trend in the data between candidate effector
number and total number of mimps identified.

Note: the correlation between variables X and Y is equal to the
correlation between variables Y and X so the order of the variables in
the `cor.test()` function does not matter.

``` r
#Correlation test:
effectors_v_mimps <- cor.test(stats_data$candidate_effectors, stats_data$mimps, )
effectors_v_mimps
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  stats_data$candidate_effectors and stats_data$mimps
    ## t = 11.712, df = 40, p-value = 1.674e-14
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.7862456 0.9340469
    ## sample estimates:
    ##       cor 
    ## 0.8799097

For the correlation between “candidate_effectors” and “mimps,” the
Pearson’s correlation coefficient is 0.88. This indicates a strong
positive correlation between these two variables. The p-value
(1.674e-14) is very small, suggesting that this correlation is
statistically significant. The 95 percent confidence interval for the
correlation coefficient ranges from 0.786 to 0.934. Overall, it implies
a robust and positive association between “candidate_effectors” and
“mimps” in my data.

``` r
#Correlation test:
effectors_v_assembly_size <- cor.test(stats_data$candidate_effectors, stats_data$assembly_size) 
effectors_v_assembly_size
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  stats_data$candidate_effectors and stats_data$assembly_size
    ## t = 13.854, df = 40, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8372508 0.9507566
    ## sample estimates:
    ##      cor 
    ## 0.909695

Similarly, for the correlation between “candidate_effectors” and
“assembly_size,” the Pearson’s correlation coefficient is 0.91. This
indicates a very strong positive correlation between these two
variables. The p-value is extremely small (\< 2.2e-16), indicating that
this correlation is highly statistically significant. The 95 percent
confidence interval for the correlation coefficient ranges from 0.837 to
0.951. Overall, it suggests a strong and positive association between
“candidate_effectors” and “assembly_size” in my dataset.

``` r
#Correlation test:
mimps_v_assembly_size <- cor.test(stats_data$mimps, stats_data$assembly_size)
mimps_v_assembly_size
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  stats_data$mimps and stats_data$assembly_size
    ## t = 11.177, df = 40, p-value = 7.091e-14
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.7701010 0.9286194
    ## sample estimates:
    ##       cor 
    ## 0.8703258

The Pearson’s correlation coefficient between “mimps” and
“assembly_size” variables is 0.87. This indicates a strong positive
correlation between the two variables. The p-value is very small
(7.091e-14), suggesting that this correlation is statistically
significant. The 95 percent confidence interval for the correlation
coefficient ranges from 0.770 to 0.929. Overall, it implies a robust and
positive association between the variables in my data.

### Relationship between predicted effectors, *Mimps* and Assembly size

Although I have ploted and calculated correlation between the number of
*mimps*, candidate, effectors and assembly size, it is not easy to
visualise/see in one place. Therefore, I generated a single plot which
contains the number of *mimps*, candidate effectors, and assembly size
for each isolate. This makes it quick and easy to visualise, and see if
there is a relationship between the variables.

First, I have to prepare my data set. I used the metadata file loaded
initially.

``` r
# ---- Prepare Data for plotting --- #
#Extract the isolate, assembly size, total number of mimps and effectors columns.
stats_plot_data <- select(metadata,"species", "species_group", "fsp", "isolate_code","genome_size","no._mimps","no._cand_effs") %>%
#Rename the columns to reduce the long titles
  dplyr::rename(isolate=isolate_code,                                    
         mimps=no._mimps,
         candidate_effectors=no._cand_effs,
         assembly_size =genome_size) %>%  
#We need to drop rows which were not included in the Maei analysis
  drop_na(candidate_effectors) %>%
#Merge/group the fsp and isolate code columns so that both can be plotted.
   unite(ID, c(fsp, isolate), sep = " ", remove = T, na.rm = T) %>%
#Merge/group the mimps and effector columns so that both can be plotted per strain/isolate.
  pivot_longer(cols = c(mimps,candidate_effectors), names_to="Legend", values_to="mimps_and_candidate_effectors") %>%  
  mutate(Legend = factor(Legend, levels=c('mimps','candidate_effectors')))
```

Once the data was prepared, I produced a plot using ggplot2.

``` r
#Build plot
#----------

#Generate scale for Assembly size data 
scale_right <- 100 / max(stats_plot_data$mimps_and_candidate_effectors)

#Build Plot
stats_plot <- ggplot(aes(x=reorder(ID, mimps_and_candidate_effectors)),   #Create X axis, which contains all strains/isolates assessed ordered by the total number of mimps and effectors 
               data = stats_plot_data)+  
  geom_bar(aes(y=mimps_and_candidate_effectors,                  #Plot the total number of mimps and effectors 
               fill = Legend),
           colour="grey20",
           position= 'dodge',                          #Ensure the bars are not stacked. 
           stat='identity')+                           #Add the mimp or predicted effector content. 
  scale_fill_manual("Legend", values=c("candidate_effectors" = "darkolivegreen", "mimps" = "#DDE0DA"), label=c("mimps", "candidate effectors"))+ #for some reason the labels have to be written the other way round...?
  facet_grid(~species_group,
             scales = "free_x",                        # Let the x axis vary across facets.
             space = "free_x",                         # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+                            # Move the facet labels to the bottom.
  geom_point(aes(y=assembly_size /scale_right,         #Plot assemble size over the top of the bar chart.
                 colour = "Assembly Size",             #Add assembly size to the legend. 
                 group = 1), 
             size = 3)+ 
  scale_colour_manual(" ", values=c("Assembly Size" = "grey50"))+
  theme_bw()+
  theme(legend.box="verticle",
        legend.title = element_blank(),
        legend.position = "bottom")+
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 12))+
  theme(axis.text = element_text(size = 12))+
  scale_y_continuous(name= "Total number of mimps \nand candidate effectors", #Increase ticks on Y axis.
                     breaks = scales::pretty_breaks(n = 10),                     #Create regular breaks in the Y axis. 
                     sec.axis = sec_axis( trans=~.*scale_right,                   #Insert second Y axis for the assembly size. Calculated using the right-hand scale.   
                                          name="Size of Assembly (Mb)", 
                                          breaks = scales::pretty_breaks(n = 15)
                     ))
# print the plot
plot(stats_plot)
```

<figure>
<img
src="AnalysisCandidateEffectorSets_files/figure-gfm/stats%20final%20plot-1.png"
alt="Figure 4:Relationship between candidate effectors, mimps and assembly size." />
<figcaption aria-hidden="true">Figure 4:Relationship between candidate
effectors, <em>mimps</em> and assembly size.</figcaption>
</figure>

``` r
#Prepare png of file 
ggsave("StatsOverview.png", width = 10, height = 10)
```

### Additional statsistics

I considered performing t-tests to see if there was a significant
difference in the number of *mimps* and candidate effectors identified
between race in each fsp, but as the number of fsp for *Fo.* fsp.
*cubense* was not even, and one of the three Race 1 isolates reported
might not be pathogenic (Foc1 60), I decided not to test for
significance - even using non-parametric tests. The sample number is too
low.

I wondered if there was a reduced number of effectors in the Foc
isolates compared to the other fsp. They do consistently seem to have
fewer identified. Further, I wondered if there was a difference in the
number of *mimps* and candidate effectors between different races of the
same fsp, even if it was not signifcant.

I think it is interesting, though we cannot do any really significant
stats on this, to look at candidate effector and mimp content among the
race groups.

I was expecting to see a fairly similar number of mimps and candidate
effectors among the TR4 isolates, as they are from a monophyletic clade
in tef and rbp2 phylos, but number of candidates varies. I am thinking
it may be due to assembly quality too? Looking at the Foa isolates is a
good way to do this, as the R4 and R2 assemblies are from the same
isolates but diff versions prepared by diff people using different
methods.

I plotted the mimp and effector content from the metadata and visualised
it using `ggplot.`

``` r
# ---- Prepare Data for plotting --- #
#Extract the isolate, assembly size, total number of mimps and effectors coloumns.
race_plot_data <- select(metadata,"species", "species_group", "fsp", "race" ,"isolate_code","genome_size","no._mimps","no._cand_effs") %>%
#Rename the columns to reduce the long titles
  dplyr::rename(isolate=isolate_code,                                    
         mimps=no._mimps,
         candidate_effectors=no._cand_effs,
         assembly_size =genome_size) %>%  
#We need to drop rows which were not included in the Maei analysis
  drop_na(candidate_effectors) %>%
#Merge/group the mimps and effector columns so that both can be plotted per strain/isolate.
  pivot_longer(cols = c(mimps,candidate_effectors), names_to="Legend", values_to="mimps_and_candidate_effectors") %>%  
  mutate(Legend = factor(Legend, levels=c('mimps','candidate_effectors')))


# ---- subset metadata ---- #
# extract only the fsp we are interested in and drop isolates which don't have a race classification
race_plot_data_subset <- subset(race_plot_data, grepl("lactucae|apii|cubense", fsp)) %>%
  drop_na(race)

# ---- funky cheats to lable the facet plot nicely --- #

#set new labels
new_labels_y <- c("mimps" = "mimps", "candidate_effectors" = "candidate effectors")
new_labels_x <- c("apii" = "Fo. fsp. apii", "cubense" = "Fo. fsp. cubense", "lactucae" = "Fo. fsp. lactucae")


# ---- build plots ---- #
theme_set(theme_pubr()) #set the ggpubr theme

# build plots for candidate effectors and mimps 
mimpsandcandeffs <- ggplot(race_plot_data_subset, aes(x=race, y=mimps_and_candidate_effectors)) + # plot race and mimp/candidate effector count
  geom_boxplot(aes(fill = Legend)) +
  facet_wrap(~ fsp, labeller = labeller(Legend = new_labels_y, fsp = new_labels_x), scales = "free_x") +  # split the plot by fsp and then mimps/candidate effectors
  labs(x = "Race", 
       y = "Total number of mimps \nand candidate effectors") +
  scale_fill_manual(values = c("candidate_effectors" = "darkolivegreen", "mimps" = "#DDE0DA"), label=c("mimps", "candidate effectors")) +
  theme_bw() +
  theme(strip.text.y = element_blank(),  # remove the side names as we have this shown in colour now. 
        panel.grid.major = element_blank(), # put the lines back in 
        legend.position = "bottom", 
        legend.title = element_blank())

#plot it 
plot(mimpsandcandeffs)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/effector%20and%20mimp%20race%20plots-1.png)<!-- -->

``` r
#save the plot
ggsave("MimpsAndCandEffs_FspOfInterest.png", width = 10, height = 10)
```

I want to put some values to this to dicuss in the text of my thesis.

``` r
# ---- subset the fsp of interest ---- #
# cubense
foc_stats_sum <- foc_stats %>%
  group_by(race) %>%
  summarize(count = n_distinct(label),
            mean_mimp = mean(no._mimps),
            min_mimp = min(no._mimps),
            max_mimp= max(no._mimps),
            mean_eff = mean(no._cand_effs),
            min_eff = min(no._cand_effs),
            max_eff = max(no._cand_effs),
            sum_eff = sum(no._cand_effs))
# lactucae
fola_stats_sum <- fola_stats %>%
  group_by(race) %>%
  summarize(count = n_distinct(label),
            mean_mimp = mean(no._mimps),
            min_mimp = min(no._mimps),
            max_mimp= max(no._mimps),
            mean_eff = mean(no._cand_effs),
            min_eff = min(no._cand_effs),
            max_eff = max(no._cand_effs),
            sum_eff = sum(no._cand_effs))
# apii 
foa_stats_sum <- foa_stats  %>%
  group_by(race) %>%
  summarize(count = n_distinct(label),
            mean_mimp = mean(no._mimps),
            min_mimp = min(no._mimps),
            max_mimp= max(no._mimps),
            mean_eff = mean(no._cand_effs),
            min_eff = min(no._cand_effs),
            max_eff = max(no._cand_effs),
            sum_eff = sum(no._cand_effs))

# print it nicely
knitr::kable(foc_stats_sum)
```

| race            | count | mean_mimp | min_mimp | max_mimp | mean_eff | min_eff | max_eff | sum_eff |
|:----------------|------:|----------:|---------:|---------:|---------:|--------:|--------:|--------:|
| Race 1          |     3 |     121.0 |       39 |      183 | 67.66667 |      45 |      95 |     203 |
| Tropical Race 4 |     6 |     134.0 |       87 |      165 | 82.83333 |      60 |     127 |     497 |
| NA              |     5 |     151.4 |      137 |      179 | 71.20000 |      45 |     107 |     356 |

``` r
knitr::kable(fola_stats_sum)
```

| race   | count | mean_mimp | min_mimp | max_mimp | mean_eff | min_eff | max_eff | sum_eff |
|:-------|------:|----------:|---------:|---------:|---------:|--------:|--------:|--------:|
| Race 1 |     3 |  546.0000 |      533 |      569 | 299.6667 |     260 |     343 |     899 |
| Race 4 |     3 |  583.6667 |      522 |      615 | 506.6667 |     482 |     548 |    1520 |

``` r
knitr::kable(foa_stats_sum)
```

| race   | count | mean_mimp | min_mimp | max_mimp | mean_eff | min_eff | max_eff | sum_eff |
|:-------|------:|----------:|---------:|---------:|---------:|--------:|--------:|--------:|
| Race 2 |     2 |       431 |      420 |      442 |    393.5 |     388 |     399 |     787 |
| Race 3 |     1 |       210 |      210 |      210 |    328.0 |     328 |     328 |     328 |
| Race 4 |     2 |       378 |      217 |      539 |    344.5 |     332 |     357 |     689 |

## Candidate effector distribution

Now I have an understanding of the effector distribution and the
phylogeny, I combined the phylogenies and effector profiles to generate
a binary presence/absence heatmap.

### Summary statistics of candidate effector clusters

As we are looking at clustered (0.65% ID, cd-hit (v…)) and filtered
sequences (SignalP (v5.06) and EffectorP (v2.0.1) extracted from BLAST
hits, the total number of candidate effectors per cluster per isolate
varies, but I want to just look at Presence/Absence. In order to do
this, I converted the heatmap data matrix to a binary data frame.
However, I don’t just want to work with the binary data, so I’ll make a
second df I can use later.

``` r
# ---- prep the binary heatmap data ---- #
# remove the row names temporarily 
rownames_mat<- data[,1]
mat_data<- as.matrix(data[,-1])
#make data frame binary 
binary_matrix <- as.matrix(mat_data)
binary_matrix[binary_matrix > 0] <- 1
#put rownames back
rownames(binary_matrix)<-rownames_mat

# ---- prep the heatmap data ---- #
# make the CEC data a matric
CEC_matrix<- as.matrix(data)
# make it a df 
CEC_df <- as.data.frame(CEC_matrix)
```

First, I wanted to look at the distribution of the these clusters in
numerical terms. How many CECs are there in each of the Fusarium genomes
I included. To do this, I counted the total number of 1s in each columns
from the bimary matrix, where each column represents a CEC. I then
merged that data with the metadata df, to keep track going forward.

``` r
# ---- summarise the distribution of candidate effector clusters ---- #
# convert the matrix to a data frame
CEC_binary_df <- as.data.frame(binary_matrix)

# count the number of columns (minus the label column) to get the total number of CECs
total_CECs <- ncol(CEC_binary_df)

# total number of candidate effector clusters per assembly (as a dataframe)
cluster_distib <- enframe(colSums(t(CEC_binary_df[-1])))
#rename the value column 
cluster_distib <- dplyr::rename(cluster_distib, no._CECs = value)
# merge the two data frames
metadata <- merge(metadata, cluster_distib, by.x = "label", by.y = "name", all.x = TRUE)

#subset the metadata ft 
CEC_metadata <- select(metadata,"species", "species_group", "fsp", "race", "isolate_code","genome_size","no._mimps","no._cand_effs", "no._CECs")

# print it nicely
knitr::kable(CEC_metadata)
```

| species         | species_group | fsp          | race            | isolate_code | genome_size | no.\_mimps | no.\_cand_effs | no.\_CECs |
|:----------------|:--------------|:-------------|:----------------|:-------------|------------:|-----------:|---------------:|----------:|
| F. graminearum  | Other         | NA           | NA              | PH-1         |        38.0 |          1 |              6 |         5 |
| F. mindanaoense | Other         | NA           | NA              | PD20-05      |          NA |         NA |             NA |        NA |
| F. sacchari     | Other         | NA           | NA              | FS66         |        47.5 |          9 |             12 |        12 |
| F. sacchari     | Other         | NA           | NA              | NRRL_66326   |        42.8 |         10 |             15 |        15 |
| Fusarium        | Other         | NA           | NA              | S16          |          NA |         NA |             NA |        NA |
| Fusarium        | Other         | NA           | NA              | S32          |          NA |         NA |             NA |        NA |
| Fusarium        | Other         | NA           | NA              | S6           |          NA |         NA |             NA |        NA |
| Fusarium        | Other         | NA           | NA              | SY-2         |        44.2 |          3 |             12 |        12 |
| F. oxysporum    | F. oxysporum  | endophyte    | np              | Fo47         |        50.4 |         70 |             81 |        53 |
| F. oxysporum    | F. oxysporum  | from rocket  | NA              | AJ174        |        62.6 |        419 |            169 |        78 |
| F. oxysporum    | F. oxysporum  | apii         | Race 2          | 207A         |        64.6 |        420 |            388 |       103 |
| F. oxysporum    | F. oxysporum  | apii         | Race 4          | 274AC        |        67.3 |        217 |            357 |        90 |
| F. oxysporum    | F. oxysporum  | apii         | Race 4          | AJ498        |        64.6 |        539 |            332 |        90 |
| F. oxysporum    | F. oxysporum  | apii         | Race 2          | AJ720        |        64.7 |        442 |            399 |       101 |
| F. oxysporum    | F. oxysporum  | apii         | Race 3          | NRRL38295    |        65.3 |        210 |            328 |        87 |
| F. oxysporum    | F. oxysporum  | cepae        | Race 2          | FoC_Fus2     |        53.4 |        325 |            359 |        84 |
| F. oxysporum    | F. oxysporum  | conglutinans | NA              | Fo5176       |        68.0 |        442 |            385 |        87 |
| F. oxysporum    | F. oxysporum  | coriandrii   | NA              | 3-2          |        65.4 |        478 |            315 |        99 |
| F. oxysporum    | F. oxysporum  | coriandrii   | NA              | AJ615        |        69.3 |        675 |            603 |       106 |
| F. oxysporum    | F. oxysporum  | coriandrii   | NA              | GL306        |          NA |         NA |             NA |        NA |
| F. oxysporum    | F. oxysporum  | cubense      | Race 1          | 160527       |        51.1 |        183 |             95 |        49 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | 58           |        48.2 |        165 |            108 |        45 |
| F. oxysporum    | F. oxysporum  | cubense      | Race 1          | 60           |        48.6 |         39 |             45 |        34 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | B2           |        52.9 |        142 |             60 |        42 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | C1HIR_9889   |        46.7 |        145 |             70 |        43 |
| F. oxysporum    | F. oxysporum  | cubense      | Race 1          | N2           |        47.7 |        141 |             63 |        42 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | NRRL_54006   |        46.6 |        105 |             70 |        43 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | Pers4        |        46.4 |         87 |             62 |        45 |
| F. oxysporum    | F. oxysporum  | cubense      | Tropical Race 4 | UK0001       |        48.6 |        160 |            127 |        46 |
| F. oxysporum    | F. oxysporum  | cubense      | NA              | VPRI44079    |        49.5 |        149 |            107 |        53 |
| F. oxysporum    | F. oxysporum  | cubense      | NA              | VPRI44081    |        47.2 |        137 |             45 |        40 |
| F. oxysporum    | F. oxysporum  | cubense      | NA              | VPRI44082    |        46.3 |        146 |             50 |        47 |
| F. oxysporum    | F. oxysporum  | cubense      | NA              | VPRI44083    |        46.3 |        146 |             50 |        47 |
| F. oxysporum    | F. oxysporum  | cubense      | NA              | VPRI44084    |        49.5 |        179 |            104 |        56 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 4          | AJ516        |        68.8 |        522 |            548 |        90 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 1          | AJ520        |        62.2 |        533 |            343 |        86 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 4          | AJ592        |        66.0 |        615 |            482 |        85 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 4          | AJ705        |        66.2 |        614 |            490 |        83 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 1          | AJ718        |        62.1 |        536 |            260 |        76 |
| F. oxysporum    | F. oxysporum  | lactucae     | Race 1          | AJ865        |        62.7 |        569 |            296 |        76 |
| F. oxysporum    | F. oxysporum  | lini         | NA              | 39_C0058     |        59.2 |        263 |            237 |        87 |
| F. oxysporum    | F. oxysporum  | lycopersici  | Race 2          | 4287         |        56.2 |        332 |            337 |        80 |
| F. oxysporum    | F. oxysporum  | matthiolae   | NA              | AJ260        |        60.3 |        301 |            209 |        76 |
| F. oxysporum    | F. oxysporum  | narcissus    | NA              | FON63        |        60.0 |        555 |            251 |        91 |
| F. oxysporum    | F. oxysporum  | niveum       | np              | 110407-3-1-1 |        49.7 |         52 |             39 |        36 |
| F. oxysporum    | F. oxysporum  | rapae        | NA              | Tf1208       |        59.8 |        377 |            278 |        82 |
| F. oxysporum    | F. oxysporum  | vasinfectum  | Race 1          | TF1          |        50.0 |        199 |             91 |        58 |

``` r
# # summarise it for chapter text
# summary(CEC_metadata)
```

Next I wanted to look at the size of the CECs, how many sequences are
there in each CEC? Do any fsp have unique CECs which are expanded?

``` r
# ---- Prep the heatmap data ---- #
# join the metadata in to make filtering and manipulation later easier
 CEC_df_all <- left_join(metadata, CEC_df, by = c("label" = "Isolate"))

# ---- group by fsp ---- #

# calculate the total number of CEs in a given CEC per fsp.
CEC_sizes <- CEC_df_all %>% 
  group_by(full_name) %>%
  mutate(across(starts_with("Cluster"),
              ~ as.numeric(as.character(.)))) %>% 
  summarize(across(starts_with("Cluster"), sum, na.rm = TRUE))

# use the CEC_sizes output to find the range for each fsp (highest number of CEs in a CEC and lowest number of CEs in a CEC, not counting 0)
CEC_ranges <- CEC_sizes %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "column") %>%
  filter(value != 0) %>%
  group_by(full_name) %>%
  summarise(num_highest = sum(value == max(value)), #  report the number of clusters with the max value, just incase it is more than one.
            highest = if_else(num_highest == 1, first(column[which.max(value)]), NA_character_), #  report the highest cluster unless it is > 1. 
            highest_value = max(value),
            num_lowest = sum(value == min(value)), #  report the number of clusters with the min value, becuase it is likely to be more than one.
            lowest = if_else(num_lowest == 1, first(column[which.min(value)]), NA_character_), #  report the lowest cluster unless it is > 1. 
            lowest_value = min(value))

# print it nicely
knitr::kable(CEC_ranges)
```

| full_name                      | num_highest | highest     | highest_value | num_lowest | lowest | lowest_value |
|:-------------------------------|------------:|:------------|--------------:|-----------:|:-------|-------------:|
| F. graminearum                 |           6 | NA          |             1 |          6 | NA     |            1 |
| F. oxysporum fsp. apii         |           1 | Cluster.219 |           195 |          3 | NA     |            1 |
| F. oxysporum fsp. cepae        |           1 | Cluster.294 |            52 |         57 | NA     |            1 |
| F. oxysporum fsp. conglutinans |           1 | Cluster.219 |            48 |         39 | NA     |            1 |
| F. oxysporum fsp. coriandrii   |           1 | Cluster.234 |           215 |         56 | NA     |            1 |
| F. oxysporum fsp. cubense      |           1 | Cluster.286 |           195 |         10 | NA     |            1 |
| F. oxysporum fsp. endophyte    |           1 | Cluster.239 |             9 |         41 | NA     |            1 |
| F. oxysporum fsp. from rocket  |           1 | Cluster.120 |            13 |         45 | NA     |            1 |
| F. oxysporum fsp. lactucae     |           1 | Cluster.234 |           628 |          8 | NA     |            1 |
| F. oxysporum fsp. lini         |           1 | Cluster.73  |            30 |         57 | NA     |            1 |
| F. oxysporum fsp. lycopersici  |           1 | Cluster.294 |            98 |         51 | NA     |            1 |
| F. oxysporum fsp. matthiolae   |           1 | Cluster.239 |            24 |         49 | NA     |            1 |
| F. oxysporum fsp. narcissus    |           1 | Cluster.239 |            23 |         55 | NA     |            1 |
| F. oxysporum fsp. niveum       |           1 | Cluster.20  |             3 |         34 | NA     |            1 |
| F. oxysporum fsp. rapae        |           1 | Cluster.73  |            33 |         45 | NA     |            1 |
| F. oxysporum fsp. vasinfectum  |           3 | NA          |             8 |         50 | NA     |            1 |
| F. sacchari                    |           9 | NA          |             2 |          9 | NA     |            1 |
| Fusarium                       |          12 | NA          |             1 |         12 | NA     |            1 |

Cluster234 is very large, with 268 CEs in lactucae - is this likley?

I wonder if Cluster234 is generally large, I’ll take a look at the
clusters to find out which are largest and smallest overall.

``` r
# now calculate the largest clusters
CEC_overall_sizes <- CEC_sizes %>%
  pivot_longer(cols = -full_name, names_to = "column") %>%
  filter(value != 0) %>%
  group_by(column) %>%
  summarise(total_value = sum(value),
            num_values = sum(value != 0)) %>%
  ungroup() %>%
  summarise(largest_column = column[which.max(total_value)],
            largest_total = max(total_value),
            num_largest = sum(total_value == max(total_value)),
            smallest_total = min(total_value),
            num_smallest = sum(total_value == min(total_value)))
            
# print it nicely
knitr::kable(CEC_overall_sizes)
```

| largest_column | largest_total | num_largest | smallest_total | num_smallest |
|:---------------|--------------:|------------:|---------------:|-------------:|
| Cluster.234    |           895 |           1 |              1 |           80 |

Cluster 234 is the largest overall, helped, no doubt, by lactucae. 895
sequences in that cluster seems considerable! I’m doubtful that’s
actually an effector - its possibly an ABC transporter which has not
been filtered out by EffectorP. It tends to do that.

There are 80 CECs with just 1 sequence, some may be genuine, some may
not.

I want to look at the distribution of the count of CEs in the CECs.

``` r
# ---- prepare the data ---- #

CEC_count <- CEC_df_all %>% 
  mutate(across(starts_with("Cluster"),
              ~ as.numeric(as.character(.)))) %>% 
  summarize(across(starts_with("Cluster"), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Cluster", values_to = "Count") # name the columns for ggplot

# ---- quick summary stats ---- #

# Calculate mean and standard deviation for normal distribution
CEC_mean_value <- mean(CEC_count$Count)
CEC_sd_value <- sd(CEC_count$Count)
CEC_range_value <- range(CEC_count$Count)

# ---- plot the CEC counts ---- #
CEC_count_plot <- ggplot(CEC_count, aes(x = Count)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "lightblue", colour = "black") +
  stat_function(fun = dnorm, args = list(mean = CEC_mean_value, sd = CEC_sd_value), color = "red", linewidth = 0.5) +
  labs(title = "Count of Clusters", x = "Count", y = "Density")

#plot it
plot(CEC_count_plot)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/plot%20of%20CE%20distrib%20in%20CEC-1.png)<!-- -->

The majority of the clusters have 0 or 1 CE in them, with only a handful
having \> 500 CEs.

### Build heatmap of CECs

Next, I build the CEC heatmap. I did originally do this alongside the
TEF phlyo, but it looks messy and unclear, plus you can see the fsp
clusters more clearly using ComplexHeatmap::pheatmap.

First, I clustered the data in the binary data matrix, so that it will
be ordered when I visualise the candidate effector heatmap.

``` r
# ---- Cluster the heatmap data ---- #
# normalisiation is mandatory for clustering, but as my data is binary - i did not normalise. 
# Compute hierarchical clustering of columns
heatmap_dat <- cluster_matrix(binary_matrix, method ="ward.D2")
```

Next, I built the plot. I started using `pheatmap` but that didn’t have
all the functionality I wanted, so I swapped to `ComplexHeatmap` with
the `pheatmap` plug in, so I didn’t have to start the heatmap again!

``` r
# ---- Add metadata to heatmap ---- #

# first, select the data we want from metadata and make the "label column in metadata the row names"
heatmap_metadf <- select(metadata,"label","no._mimps","no._cand_effs", "no._CECs", "race", "full_name", "isolate_code") %>%
  dplyr::rename(name=label,
                Species =full_name,
                CECs=no._CECs,
                CEs = no._cand_effs,
                mimps =no._mimps,) %>% 
  drop_na(CECs) %>%
  mutate(across('Species', str_replace, 'F. oxysporum', 'Fo.')) %>% # shorten the species name for Fusarium oxysporum.
  unite(ID, c(Species, isolate_code), sep = " ", remove = F, na.rm = T) %>%
  remove_rownames %>% 
  tibble::column_to_rownames(var="name")
# because complex heatmap cant cope with the order of the df being different from the matrix, I have to reorder out df to match the matrix.
heatmap_metadf_ordered<- heatmap_metadf[rownames(binary_matrix), ]

#add isolate ids
ID <- as.list(heatmap_metadf_ordered$ID)

# ---- set my colours ---- #

anno_colours_r <- list(
           Species = c(
             "Fo. fsp. apii" = "blue",
             "Fo. fsp. cepae" = "purple",
             "Fo. fsp. conglutinans" = "goldenrod",
             "Fo. fsp. coriandrii" = "grey90",
             "Fo. fsp. cubense" = "gold",
             "Fo. fsp. lactucae" = "darkolivegreen3",
             "Fo. fsp. lini" = "indianred",
             "Fo. fsp. lycopersici" = "tomato",
             "Fo. fsp. matthiolae" = "lavender",
             "Fo. fsp. narcissus" = "yellow",
             "Fo. fsp. niveum" = "palegreen4",
             "Fo. fsp. rapae" = "slateblue",
             "Fo. fsp. from rocket" = "tan",
             "Fo. fsp. vasinfectum" = "steelblue",
             "Fo. fsp. endophyte" = "brown",
             "F. graminearum" = "burlywood4",
             "F. sacchari" = "pink",
             "Fusarium" = "red"),
           Race = c(
             "np" = "bisque",
             "Race 1" = "gold4",
             "Race 2" = "yellow3",
             "Race 3" = "lightblue",
             "Race 4" = "navy",
             "Tropical Race 4" = "green4"))

anno_colours_cec = colorRamp2(c(0, 150), c("white", "darkorchid4"))
anno_colours_ce = colorRamp2(c(0, 600), c("white", "darkolivegreen"))
anno_colours_mimp = colorRamp2(c(0, 800), c("white", "goldenrod1"))

# ---- create heatmap annotations ---- #

# add column barplot
column_anno = HeatmapAnnotation(
  "CEC size" = anno_barplot(colSums(binary_matrix), 
                            outline = FALSE, 
                            gp = gpar(fill = "black")))

#add row data - CEC, CE and mimp count.
row_anno_l <- rowAnnotation(
  "Total CECs" = heatmap_metadf_ordered$CECs, 
  "Total CEs" = heatmap_metadf_ordered$CEs,
  "Total mimps" = heatmap_metadf_ordered$mimps,
  col = list("Total CECs" = anno_colours_cec,
             "Total CEs" = anno_colours_ce,
             "Total mimps" = anno_colours_mimp),
  gp = gpar(col = "white")
)


row_anno_r = rowAnnotation(
  "Species" = heatmap_metadf_ordered$Species, 
  "Race" = heatmap_metadf_ordered$race, 
  col = anno_colours_r, 
  na_col = NA, 
  gp = gpar(col = "white"))
 
# ---- pheatmap plot of CECs ---- #

effector_heatmap <- ComplexHeatmap::pheatmap(
  binary_matrix, 
  color = colorRampPalette(c("grey90", "black"))(2), 
  name = "Binary distribution",
  legend = T,
  heatmap_legend_param = list(
    #at = seq(1, 10, by = 1),  #wär gleich: at = 1:10,       
    at = 0:1,
    legend_gp = gpar(fill = 0:1, fontsize = 2),
    color_bar = "discrete"
  ),
  legend_labels = c("Absent", "Present"),
  show_colnames = F,
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  #cellwidth = 2, 
  #cellheight = 20,
  border_color = NA,
  treeheight_row = 80,
  treeheight_col = 20,
  na_col = "white", 
  row_labels = ID,
  row_names_side = "right",
  top_annotation = column_anno,
  left_annotation = row_anno_l,
  right_annotation = row_anno_r,
  fontsize = 12
)

# add bar plot annotations
effector_heatmap 
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/final%20heatmap-1.png)<!-- -->

``` r
#save it
png(file="EffectorsHeatmap.png", width = 24, height = 12, unit = "in", res = 150)
draw(effector_heatmap)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

As well as plotting the CEC distribution, I also want to get a numerical
overview. First, I looked at the overall distribution and summary
statistics of CECs.

``` r
# ---- candidate effector clusters in Fo ---- #

# extract just the Fo rows from the CEC_metadata
Fo_CEC_medtadata <- subset(CEC_metadata, !grepl("Other", species_group))
#summarise
summary(Fo_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:39          Length:39          Length:39          Length:39         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size      no._mimps   no._cand_effs   
    ##  Length:39          Min.   :46.30   Min.   : 39   Min.   : 39.00  
    ##  Class :character   1st Qu.:48.83   1st Qu.:146   1st Qu.: 72.75  
    ##  Mode  :character   Median :57.70   Median :240   Median :223.00  
    ##                     Mean   :56.69   Mean   :306   Mean   :227.18  
    ##                     3rd Qu.:64.60   3rd Qu.:469   3rd Qu.:341.50  
    ##                     Max.   :69.30   Max.   :675   Max.   :603.00  
    ##                     NA's   :1       NA's   :1     NA's   :1       
    ##     no._CECs     
    ##  Min.   : 34.00  
    ##  1st Qu.: 46.25  
    ##  Median : 76.00  
    ##  Mean   : 68.84  
    ##  3rd Qu.: 87.00  
    ##  Max.   :106.00  
    ##  NA's   :1

The number of CECs in Fo. only doesn’t vary as widely as the number of
CECs across all Fusarium assemblies included (see range).

Next, I looked at specific fsp. of interest as well as comparing the *F.
sacchari* genome assemblies.

``` r
# ---- candidate effector clusters in fsp ---- #

# extract just the cubense rows from the CEC_metadata
Foc_CEC_medtadata <- subset(CEC_metadata, grepl("cubense", fsp))
summary(Foc_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:14          Length:14          Length:14          Length:14         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size      no._mimps     no._cand_effs   
    ##  Length:14          Min.   :46.30   Min.   : 39.0   Min.   : 45.00  
    ##  Class :character   1st Qu.:46.62   1st Qu.:138.0   1st Qu.: 52.50  
    ##  Mode  :character   Median :47.95   Median :145.5   Median : 66.50  
    ##                     Mean   :48.26   Mean   :137.4   Mean   : 75.43  
    ##                     3rd Qu.:49.27   3rd Qu.:157.2   3rd Qu.:101.75  
    ##                     Max.   :52.90   Max.   :183.0   Max.   :127.00  
    ##     no._CECs    
    ##  Min.   :34.00  
    ##  1st Qu.:42.25  
    ##  Median :45.00  
    ##  Mean   :45.14  
    ##  3rd Qu.:47.00  
    ##  Max.   :56.00

``` r
# extract just the lactucae rows from the CEC_metadata
Fola_CEC_medtadata <- subset(CEC_metadata, grepl("lactucae", fsp))
summary(Fola_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:6           Length:6           Length:6           Length:6          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size      no._mimps     no._cand_effs  
    ##  Length:6           Min.   :62.10   Min.   :522.0   Min.   :260.0  
    ##  Class :character   1st Qu.:62.33   1st Qu.:533.8   1st Qu.:307.8  
    ##  Mode  :character   Median :64.35   Median :552.5   Median :412.5  
    ##                     Mean   :64.67   Mean   :564.8   Mean   :403.2  
    ##                     3rd Qu.:66.15   3rd Qu.:602.8   3rd Qu.:488.0  
    ##                     Max.   :68.80   Max.   :615.0   Max.   :548.0  
    ##     no._CECs    
    ##  Min.   :76.00  
    ##  1st Qu.:77.75  
    ##  Median :84.00  
    ##  Mean   :82.67  
    ##  3rd Qu.:85.75  
    ##  Max.   :90.00

``` r
# extract just the apii rows from the CEC_metadata
Foa_CEC_medtadata <- subset(CEC_metadata, grepl("apii", fsp))
summary(Foa_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:5           Length:5           Length:5           Length:5          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size     no._mimps     no._cand_effs  
    ##  Length:5           Min.   :64.6   Min.   :210.0   Min.   :328.0  
    ##  Class :character   1st Qu.:64.6   1st Qu.:217.0   1st Qu.:332.0  
    ##  Mode  :character   Median :64.7   Median :420.0   Median :357.0  
    ##                     Mean   :65.3   Mean   :365.6   Mean   :360.8  
    ##                     3rd Qu.:65.3   3rd Qu.:442.0   3rd Qu.:388.0  
    ##                     Max.   :67.3   Max.   :539.0   Max.   :399.0  
    ##     no._CECs    
    ##  Min.   : 87.0  
    ##  1st Qu.: 90.0  
    ##  Median : 90.0  
    ##  Mean   : 94.2  
    ##  3rd Qu.:101.0  
    ##  Max.   :103.0

``` r
# extract just the coriandrii rows from the CEC_metadata
Foci_CEC_medtadata <- subset(CEC_metadata, grepl("coriandrii", fsp))
summary(Foci_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:3           Length:3           Length:3           Length:3          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size      no._mimps     no._cand_effs
    ##  Length:3           Min.   :65.40   Min.   :478.0   Min.   :315  
    ##  Class :character   1st Qu.:66.38   1st Qu.:527.2   1st Qu.:387  
    ##  Mode  :character   Median :67.35   Median :576.5   Median :459  
    ##                     Mean   :67.35   Mean   :576.5   Mean   :459  
    ##                     3rd Qu.:68.33   3rd Qu.:625.8   3rd Qu.:531  
    ##                     Max.   :69.30   Max.   :675.0   Max.   :603  
    ##                     NA's   :1       NA's   :1       NA's   :1    
    ##     no._CECs    
    ##  Min.   : 99.0  
    ##  1st Qu.:100.8  
    ##  Median :102.5  
    ##  Mean   :102.5  
    ##  3rd Qu.:104.2  
    ##  Max.   :106.0  
    ##  NA's   :1

``` r
# ---- candidate effector clusters in F. sacchari and SY-2 ---- #

# extract just the cubense rows from the CEC_metadata
Fs_CEC_medtadata <- subset(CEC_metadata, grepl("FS66|NRRL_66326|SY-2", isolate_code))
summary(Fs_CEC_medtadata)
```

    ##    species          species_group          fsp                race          
    ##  Length:3           Length:3           Length:3           Length:3          
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  isolate_code        genome_size      no._mimps      no._cand_effs 
    ##  Length:3           Min.   :42.80   Min.   : 3.000   Min.   :12.0  
    ##  Class :character   1st Qu.:43.50   1st Qu.: 6.000   1st Qu.:12.0  
    ##  Mode  :character   Median :44.20   Median : 9.000   Median :12.0  
    ##                     Mean   :44.83   Mean   : 7.333   Mean   :13.0  
    ##                     3rd Qu.:45.85   3rd Qu.: 9.500   3rd Qu.:13.5  
    ##                     Max.   :47.50   Max.   :10.000   Max.   :15.0  
    ##     no._CECs   
    ##  Min.   :12.0  
    ##  1st Qu.:12.0  
    ##  Median :12.0  
    ##  Mean   :13.0  
    ##  3rd Qu.:13.5  
    ##  Max.   :15.0

#### Core CEC distribution

I also wanted to know what CECs were shared among all assemblies, and
which were shared among all Fo. and which were shared among all
assemblies in a specific fsp.

``` r
# ---- summarise the distribution of candidate effector clusters in all assemblies  ---- #

# count the number of rows where the column total is >= the number of rows (a shared candidate effector cluster!)
shared_cluster <- CEC_binary_df %>%
  select_if(colSums(CEC_binary_df) >= nrow(CEC_binary_df))
# count the number of columns 
ncol(shared_cluster)
```

    ## [1] 1

``` r
# ---- summarise the distribution of candidate effector clusters in Fo ---- #

# subset only the Fo rows
heatmap_Fo_only <- subset(CEC_binary_df, grepl("^Fo", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Fo_shared <- heatmap_Fo_only %>%
  select_if(colSums(heatmap_Fo_only) >= nrow(heatmap_Fo_only))
# count the number of columns 
ncol(heatmap_Fo_shared)
```

    ## [1] 8

``` r
# ---- summarise the distribution of candidate effector clusters in cubense ---- #

# subset only the Foc rows
heatmap_Foc_only <- subset(CEC_binary_df, grepl("^Fo._fsp._cubense", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foc_shared <- heatmap_Foc_only %>%
  select_if(colSums(heatmap_Foc_only) >= nrow(heatmap_Foc_only))
# count the number of columns 
ncol(heatmap_Foc_shared)
```

    ## [1] 16

``` r
# ---- summarise the distribution of candidate effector clusters in lactucae ---- #

# subset only the Fola rows
heatmap_Fola_only <- subset(CEC_binary_df, grepl("^Fo._fsp._lactucae", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Fola_shared <- heatmap_Fola_only %>%
  select_if(colSums(heatmap_Fola_only) >= nrow(heatmap_Fola_only))
# count the number of columns 
ncol(heatmap_Fola_shared)
```

    ## [1] 56

As apii and coriandrii are have similar genomes, and apii r4 can infect
coriander, I want to examine their shared CECs more closely.

First, how many core CECs in apii, how many core CECs shared between
apii R3 and R4 and what are the CECs that are not shared?

``` r
# ---- summarise the distribution of candidate effector clusters in apii ---- #

# subset only the Foa rows
heatmap_Foa_only <- subset(CEC_binary_df, grepl("^Fo._fsp._apii", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foa_shared <- heatmap_Foa_only %>%
  select_if(colSums(heatmap_Foa_only) >= nrow(heatmap_Foa_only))
# count the number of columns 
ncol(heatmap_Foa_shared)
```

    ## [1] 54

``` r
# ---- summarise the distribution of candidate effector clusters in apii R3 and R4 (not R2)---- #

# subset only the Foa rows
heatmap_Foa_r3r4_only <- subset(CEC_binary_df, grepl("^Fo._fsp._apii_274.AC|^Fo._fsp._apii_AJ498|^Fo._fsp._apii_NRRL38295", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foa_r3r4_only_shared <- heatmap_Foa_r3r4_only %>%
  select_if(colSums(heatmap_Foa_r3r4_only) >= nrow(heatmap_Foa_r3r4_only))
# count the number of columns 
ncol(heatmap_Foa_r3r4_only_shared)
```

    ## [1] 85

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foa race 3/race4 specific matrix
heatmap_Foa_r3r4_only_drop_0 <- heatmap_Foa_r3r4_only[, colSums(heatmap_Foa_r3r4_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foa_r3r4_only_drop_0), names(heatmap_Foa_r3r4_only_shared))
```

Next, I want to look at coriandrii alone. How many core CECs between
coriandrii and how many core CECs when you include apii?

``` cec

# ---- summarise the distribution of candidate effector clusters in coriandrii ---- #

# subset only the Foci rows
heatmap_Foci_only <- subset(CEC_binary_df, grepl("^Fo._fsp._coriandrii", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foci_shared <- heatmap_Foci_only %>%
  select_if(colSums(heatmap_Foci_only) >= nrow(heatmap_Foci_only))
# count the number of columns 
ncol(heatmap_Foci_shared)

# ---- summarise the distribution of candidate effector clusters in apii and coridanrii (as they share some hosts) ---- #

# subset only the Foa and Foci rows
heatmap_Foa_c_only <- subset(CEC_binary_df, grepl("^Fo._fsp._apii|^Fo._fsp._coriandrii", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foa_c_shared <- heatmap_Foa_c_only %>%
  select_if(colSums(heatmap_Foa_c_only) >= nrow(heatmap_Foa_c_only))
# count the number of columns 
ncol(heatmap_Foa_c_shared)
```

Knowing that coriandrii 3-2 and apii race 3 and race 4 are closely
related and have a similar CEC profile, how many CECs are shared between
just these isolates

``` r
# ---- Foa monophyletic group CEC counts ---- #
# subset only the Foa rows
heatmap_Foa_r3r4_andcorr_only <- subset(CEC_binary_df, grepl("^Fo._fsp._apii_274.AC|^Fo._fsp._apii_AJ498|^Fo._fsp._apii_NRRL38295|Fo._fsp._coriandrii_3-2", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foa_r3r4_andcorr_only_shared <- heatmap_Foa_r3r4_andcorr_only %>%
  select_if(colSums(heatmap_Foa_r3r4_andcorr_only) >= nrow(heatmap_Foa_r3r4_andcorr_only))
# count the number of columns 
ncol(heatmap_Foa_r3r4_andcorr_only_shared)
```

    ## [1] 73

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foa race 3/race4 specific matrix
heatmap_Foa_r3r4_andcorr_only_drop_0 <- heatmap_Foa_r3r4_andcorr_only[, colSums(heatmap_Foa_r3r4_andcorr_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foa_r3r4_andcorr_only_drop_0), names(heatmap_Foa_r3r4_andcorr_only_shared))
```

What about apii R2 and the other coriandrii isolate?

``` r
# ---- Foa monophyletic group CEC counts ---- #
# subset only the Foa rows
heatmap_Foa_r2_and_corr_only <- subset(CEC_binary_df, grepl("^Fo._fsp._apii_207.A|^Fo._fsp._apii_AJ720|Fo._fsp._coriandrii_AJ615", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foa_r2_and_corr_only_shared <- heatmap_Foa_r2_and_corr_only %>%
  select_if(colSums(heatmap_Foa_r2_and_corr_only) >= nrow(heatmap_Foa_r2_and_corr_only))
# count the number of columns 
ncol(heatmap_Foa_r2_and_corr_only_shared)
```

    ## [1] 67

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foa race 3/race4 specific matrix
heatmap_Foa_r2_and_corr_only_drop_0 <- heatmap_Foa_r2_and_corr_only[, colSums(heatmap_Foa_r2_and_corr_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foa_r2_and_corr_only_drop_0), names(heatmap_Foa_r2_and_corr_only_shared))
```

Next, how many CECs are shared among all of the Fs and Fusarium from
TNAU.

``` r
# ---- summarise the distribution of candidate effector clusters in F. sacchari and SY-2) ---- #

# subset only the Foa and Foci rows
heatmap_Fs_only <- subset(CEC_binary_df, grepl("^F._sacchari_|^F._TNAU", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Fs_shared <- heatmap_Fs_only %>%
  select_if(colSums(heatmap_Fs_only) >= nrow(heatmap_Fs_only))
# count the number of columns 
ncol(heatmap_Fs_shared)
```

    ## [1] 9

Are the ‘core CECs’ in Fs also found in Focub?

``` r
# ---- summarise the distribution of candidate effector clusters in F. sacchari and SY-2) ---- #

# subset only the Foa and Foci rows
heatmap_banana_only <- subset(CEC_binary_df, grepl("^Fo._fsp._cubense|^F._sacchari_|^F._TNAU", rownames(CEC_binary_df)))
# count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_banana_shared <- heatmap_banana_only %>%
  select_if(colSums(heatmap_banana_only) >= nrow(heatmap_banana_only))
# count the number of columns 
ncol(heatmap_banana_shared)
```

    ## [1] 4

How are CECs distributed among races?

``` r
# ---- subset the fsp of interest ---- #
# cubense
foc_CEC_stats_sum <- Foc_CEC_medtadata %>%
  group_by(race) %>%
  summarize(count = n_distinct(isolate_code),
            mean_CEC = mean(no._CECs),
            min_CEC = min(no._CECs),
            max_CEC = max(no._CECs))
# lactucae
fola_CEC_stats_sum <- Fola_CEC_medtadata %>%
  group_by(race) %>%
  summarize(count = n_distinct(isolate_code),
            mean_CEC = mean(no._CECs),
            min_CEC = min(no._CECs),
            max_CEC= max(no._CECs))
# apii 
foa_CEC_stats_sum <- Foa_CEC_medtadata  %>%
  group_by(race) %>%
  summarize(count = n_distinct(isolate_code),
            mean_CEC = mean(no._CECs),
            min_CEC = min(no._CECs),
            max_CEC = max(no._CECs))

# print it nicely
knitr::kable(foc_CEC_stats_sum)
```

| race            | count | mean_CEC | min_CEC | max_CEC |
|:----------------|------:|---------:|--------:|--------:|
| Race 1          |     3 | 41.66667 |      34 |      49 |
| Tropical Race 4 |     6 | 44.00000 |      42 |      46 |
| NA              |     5 | 48.60000 |      40 |      56 |

``` r
knitr::kable(fola_CEC_stats_sum)
```

| race   | count | mean_CEC | min_CEC | max_CEC |
|:-------|------:|---------:|--------:|--------:|
| Race 1 |     3 | 79.33333 |      76 |      86 |
| Race 4 |     3 | 86.00000 |      83 |      90 |

``` r
knitr::kable(foa_CEC_stats_sum)
```

| race   | count | mean_CEC | min_CEC | max_CEC |
|:-------|------:|---------:|--------:|--------:|
| Race 2 |     2 |      102 |     101 |     103 |
| Race 3 |     1 |       87 |      87 |      87 |
| Race 4 |     2 |       90 |      90 |      90 |

### Candidate effector cluster distribution in Fo. fsp. cubense

Now I have my overall heatmap, I want to look at some of the fsp in more
detail, including the TEF phylogeny data - particularly those for which
we have multiple races available. I subset my metadata, searched for all
cases that did not match the regex “Fo.\_fsp.\_cubense” and dropped them
from the tree using the `drop.tip` function. I the reconstructed my tree
using just the Foc tef phylo and foc effector profiles using the
gheatmap package from ggtree.

``` r
# ---- subset foc metadata ---- #
#identify all rows in the metadata which do not contain cubense 
foc_set_df <- subset(metadata, !grepl("Fo._fsp._cubense|F._sacchari|SY-2|Fo47", label))
#subset just tip labels
foc_set <- data.frame("label" = foc_set_df[,c("label")])
# convert it to a list 
foc_set_2 <- paste(foc_set$label, sep = ",")

# subset race data
foc_set_df_race <- data.frame("race" = metadata[,c("race")] )
rownames(foc_set_df_race) <- metadata$label

# ---- subset foc heatmap data ---- #
# reduce the white space in the heatmap but filtering columns where there is no data for foc
# first we extract only the foc rows using the same approach as for the metadata, but instead we perform on the binary matrix
foc_heat_df <- subset(binary_matrix, grepl("Fo._fsp._cubense|F._sacchari|SY-2|Fo47", rownames(binary_matrix)))
# now we need to drop the empty columns 
foc_heat_df <- foc_heat_df[, colSums(foc_heat_df != 0) > 0]

node_data <- data.frame(node=c(20, 34), type=c("FOSC", "FSSC"))

# ---- Cluster the heatmap data (again)---- #
# normalisiation is mandatory for clustering, but as my data is binary - i did not normalise. 
# Compute hierarchical clustering of columns
foc_heat_df <- cluster_matrix(foc_heat_df, dim = 'col', method ="ward.D2")

# ---- build the basic foc tree ---- #
#dop the non foc tip from the tree
tree_reduced <- drop.tip(tree, foc_set_2)
#build a tree from this data and add in out metadata
foc_tree <- ggtree(tree_reduced, ladderize = T ) %<+% metadata 

# ---- build full foc specific tree ---- #
# now make the tree pretty 
foc_tree_2 <- foc_tree +  
  geom_tiplab(aes(label = isolate_code), color = "grey20", offset = 0.0018, linetype = "blank", geom = "text", align = T) +
  geom_hilight(node=20, fill="gold") +
  geom_hilight(node=34, fill="lightpink") +
  geom_tippoint(aes(shape = source), size = 3) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 80), nudge_x = -0.0027) +
  geom_rootedge() +
  theme(legend.position = "bottom")

# ---- add race data ---- #
#add extra scale so we can plot race with colour
foc_tree_3 <- foc_tree_2 + new_scale_fill()

foc_tree_4 <- gheatmap(foc_tree_3, foc_set_df_race, 
               offset = 0.015, 
               width = 0.06,
               color = "grey",
               colnames = FALSE) +
  scale_fill_manual(name = "Race",
                    values = c("bisque","gold4", "green4"), na.value = "white") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

# ---- add candidate effector clusters ---- #

# add add scale for candidate effector clusters
foc_tree_5 <- foc_tree_4 + new_scale_fill()

# plot the candidate effector clusters
foc_tree_6 <- gheatmap(foc_tree_5, foc_heat_df, offset=0.02, colnames=T, colnames_angle=90, hjust=0.6, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 4)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "grey90", high = "black",
                        breaks = c("Absent","Present"),
                        na.value = "grey")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
   geom_treescale(x=0, y=0.1, width = 0.004) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
    theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 8,
                                   colour = "white" )) +
  coord_cartesian(clip = "off")

# plot the tree
plot(foc_tree_6)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/cubense%20effector%20distib-1.png)<!-- -->

``` r
# save the tree
ggsave("HeatmapAndPhylo_BananaPathOnly.png", width = 20, height = 10)
```

We can see that the effector profiles differ based on race and overall
tef1 phylogeny. Interestingly there is one candidate which is shared
across the R1 and suspected STR4 isolates but not found in the TR4
isolates - cognate R gene?

Further, Foc1_60 displays quite a reduced CEC profile compared to the
other R1 isolates, inlcuding N2 (from the same tef lineage).

I’ll just quickly plot total CEC count per isolate to see if it as an
outlier.

``` r
# ---- prepare the data ---- #
#identify all rows in the metadata which do not contain cubense 
foc_stats_plot_data  <- subset(metadata, grepl("Fo._fsp._cubense|Fo._fsp._apii|Fo._fsp._lactucae|Fo._fsp._coriandrii_|F._sacchari|SY-2|Fo._fsp._niveum|Fo._Fo47", label)) %>%
  drop_na(no._CECs) 
#subset just tip labels

new_labels <- c("F. oxysporum fsp. coriandrii" = "Fo. fsp.\ncoriandrii","F. oxysporum fsp. apii" = "Fo. fsp.\napii", "F. oxysporum fsp. lactucae" = "Fo. fsp.\nlactucae","F. oxysporum fsp. cubense" = "Fo. fsp.\ncubense", "F. oxysporum fsp. endophyte" = "Fo.\nendo-\nphyte", "F. oxysporum fsp. niveum" = "Fo. fsp.\nniveum", "F. sacchari" = "F. sacchari" , "Fusarium" = "un-\nknown" )

# ---- build the plot ---- #
#Generate scale for Assembly size data 

foc_stats_plot <- ggplot(foc_stats_plot_data, aes(x= reorder(isolate_code, no._CECs), y=no._CECs)) + # plot race and mimp/candidate effector count
  geom_point(aes(colour = race, 
                 size = genome_size),
           position= 'dodge',                          #Ensure the bars are not stacked. 
           stat='identity')+      #Add the mimp or predicted effector content. 
  scale_color_manual("Race", values=c(
             "np" = "bisque",
             "Race 1" = "gold2",
             "Race 2" = "yellow3",
             "Race 3" = "lightblue",
             "Race 4" = "navy",
             "Tropical Race 4" = "green4"), 
             na.value = "grey20") + #for some reason the labels have to be written the other way round...?
  scale_size_binned("Genome\nassembly size") +
  facet_grid( ~ reorder(full_name, no._CECs),
              labeller = as_labeller(new_labels),  # becuase the reorder has changed the name of the faceting variable, i have to used as_labeller instead!
             scales = "free_x",                        # Let the x axis vary across facets.
             space = "free_x",                         # Let the width of facets vary and force all bars to have the same width.
             switch = "x") +                            # Move the facet labels to the bottom.
  theme_bw()+
  theme(legend.box="hortizontal",
        legend.position = "right")+
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 12))+
  theme(axis.text = element_text(size = 12))+
  theme(strip.text = element_text(size = 10)) + # shrink the text because it keeps getting cut
  scale_y_continuous(name= "Total number of Candidate\nEffector Clusters") #Increase ticks on Y axis.
                    
#plot it 
plot(foc_stats_plot)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/plot%20of%20CEC%20count-1.png)<!-- -->

``` r
# save it
ggsave("./CecDistribinFspOfInterest.png", width = 16, height = 7)
```

What does the distribution of clusters look like in terms of numbers?
How many are shared among the differnt races of Foc?

``` r
# ---- Foc tr4 monophyletic group CEC counts ---- #
# subset only the Foc tr4 rows
heatmap_Foc_TR4_only <- subset(CEC_binary_df, grepl("Fo._fsp._cubense_B2|Fo._fsp._cubense_C1HIR_9889|Fo._fsp._cubense_58|Fo._fsp._cubense_NRRL_54006|Fo._fsp._cubense_Pers4|Fo._fsp._cubense_UK0001", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foc_TR4_only_shared <- heatmap_Foc_TR4_only %>%
  select_if(colSums(heatmap_Foc_TR4_only) >= nrow(heatmap_Foc_TR4_only))
# count the number of columns 
ncol(heatmap_Foc_TR4_only_shared)
```

    ## [1] 38

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foc tr4 specific matrix
heatmap_Foc_TR4_only_drop_0 <- heatmap_Foc_TR4_only[, colSums(heatmap_Foc_TR4_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foc_TR4_only_drop_0), names(heatmap_Foc_TR4_only_shared))
```

``` r
# ---- Foc tr4 monophyletic group CEC counts ---- #
# subset only the Foc tr4 rows
heatmap_Foc_TR4_only <- subset(CEC_binary_df, grepl("Fo._fsp._cubense_B2|Fo._fsp._cubense_C1HIR_9889|Fo._fsp._cubense_58|Fo._fsp._cubense_NRRL_54006|Fo._fsp._cubense_Pers4|Fo._fsp._cubense_UK0001|Fo._fsp._cubense_VPRI44081|Fo._fsp._cubense_VPRI44082|Fo._fsp._cubense_VPRI44083", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foc_TR4_only_shared <- heatmap_Foc_TR4_only %>%
  select_if(colSums(heatmap_Foc_TR4_only) >= nrow(heatmap_Foc_TR4_only))
# count the number of columns 
ncol(heatmap_Foc_TR4_only_shared)
```

    ## [1] 27

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foc tr4 specific matrix
heatmap_Foc_TR4_only_drop_0 <- heatmap_Foc_TR4_only[, colSums(heatmap_Foc_TR4_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foc_TR4_only_drop_0), names(heatmap_Foc_TR4_only_shared))
```

``` r
# ---- Foc tr4 monophyletic group CEC counts ---- #
# subset only the Foc tr4 rows
heatmap_Foc_TR4_only <- subset(CEC_binary_df, grepl("Fo._fsp._cubense_VPRI44081|Fo._fsp._cubense_VPRI44082|Fo._fsp._cubense_VPRI44083", rownames(CEC_binary_df)))
  # count the number of rows where the column total is >= the number of Fo rows (a shared candidate effector cluster!)
heatmap_Foc_TR4_only_shared <- heatmap_Foc_TR4_only %>%
  select_if(colSums(heatmap_Foc_TR4_only) >= nrow(heatmap_Foc_TR4_only))
# count the number of columns 
ncol(heatmap_Foc_TR4_only_shared)
```

    ## [1] 38

``` r
# -- identify which CECs are unique 
# first drop all the columns where all values == 0 in the Foc tr4 specific matrix
heatmap_Foc_TR4_only_drop_0 <- heatmap_Foc_TR4_only[, colSums(heatmap_Foc_TR4_only != 0) > 0]
#next use set diff to identify the different columns 
not_shared_cols <- setdiff(names(heatmap_Foc_TR4_only_drop_0), names(heatmap_Foc_TR4_only_shared))
```

### Candidate effector cluster distribution in Fo. fsp. lactucae and Fo. fsp. apii

``` r
# ---- subset metadata ---- #
#identify all rows in the metadata which do not contain apii and coriandrii 
foa_c_set_df <- subset(metadata, !grepl("Fo._fsp._coriandrii|Fo._fsp._apii|Fo._Fo47", label))
#subset just tip labels
foa_c_set <- data.frame("label" = foa_c_set_df[,c("label")])
# convert it to a list 
foa_c_set_2 <- paste(foa_c_set$label, sep = ",")

# subset race data
foa_c_set_df_race <- data.frame("race" = metadata[,c("race")] )
rownames(foa_c_set_df_race) <- metadata$label

# ---- subset heatmap data ---- #
# reduce the white space in the heatmap but filtering columns where there is no data for foc
# first we extract only the foc rows using the same approach as for the metadata, but instead we perform on the binary matrix
foa_c_heat_df <- subset(binary_matrix, grepl("Fo._fsp._coriandrii|Fo._fsp._apii|Fo._Fo47", rownames(binary_matrix)))
# now we need to drop the empty columns 
foa_c_heat_df <- foa_c_heat_df[, colSums(foa_c_heat_df != 0) > 0]

# ---- Cluster the heatmap data (again)---- #
# normalisiation is mandatory for clustering, but as my data is binary - i did not normalise. 
# Compute hierarchical clustering of columns
foa_c_heat_df <- cluster_matrix(foa_c_heat_df, dim = 'col', method ="ward.D2")

# ---- build the basic foc tree ---- #
#dop the non foc tip from the tree
tree_reduced <- drop.tip(tree, foa_c_set_2)
#build a tree from this data and add in out metadata
foa_c_tree <- ggtree(tree_reduced, ladderize = T ) %<+% metadata 

# ---- build full foc specific tree ---- #
# now make the tree pretty 
foa_c_tree_2 <- foa_c_tree +  
  #geom_tiplab(aes(label = full_name), offset = 0.00001) +
  geom_tiplab(aes(label = fsp), color = "grey20", offset = 0.0003, linetype = "blank", geom = "text", align = F) +
  geom_tiplab(aes(label = isolate_code), color = "grey20", offset = 0.0025, linetype = "blank", geom = "text", align = TRUE) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 60), nudge_x = -0.0008) +
  geom_tippoint(aes(shape = source), size = 3) +
  geom_rootedge() +
  theme(legend.position = "bottom")
# ---- add race data ---- #
#add extra scale so we can plot race with colour
foa_c_tree_3 <- foa_c_tree_2 + new_scale_fill()

foa_c_tree_4 <- gheatmap(foa_c_tree_3, foa_c_set_df_race, 
               offset = 0.006, 
               width = 0.1,
               color = "grey",
               colnames = FALSE) +
  scale_fill_manual(name = "Race",
                    values = c("bisque", "yellow3", "lightblue" ,"navy"), na.value = "white") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

# ---- add candidate effector clusters ---- #

#add extra scale so we can plot race with colour
foa_c_tree_5 <- foa_c_tree_4 + new_scale_fill()

# add race data
foa_c_tree_6 <- gheatmap(foa_c_tree_5, foa_c_heat_df, offset=0.007, colnames=T, colnames_angle=90, hjust=0.4, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 10)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "grey90", high = "black",
                        breaks = c("Absent","Present"),
                        na.value = "white")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
   geom_treescale(x=0, y=.1, width = 0.004) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
    theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 8,
                                   colour = "white" )) +
  coord_cartesian(clip = "off") 
# plot the tree
plot(foa_c_tree_6)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/coriandrii%20and%20apii%20effector%20distib-1.png)<!-- -->

``` r
# save the tree
ggsave("HeatmapAndPhylo_ApiiAndCoriandriiOnly.png", width = 20, height = 10)
```

As I compared effector content between races for Foc, I also want to do
the same thing for Fola. Again, I would have liked to performed some
statistical analysis, but pretty much all blogs, posts, and papers I
have come across have recommended not doing it with a such a small
sample size. I have therefore just plotted the effector data for the
Fola isolates. I also included matthiolae, as it appears to be in the
same clade based on TEF, but displays a different candidate effector
profiles.

``` r
# ---- subset fola metadata ---- #
#identify all rows in the metadata which do not contain lactucae 
fola_set_df <- subset(metadata, !grepl("Fo._fsp._lactucae|Fo._fsp._matthiolae|Fo._Fo47", label))
#subset just tip labels
fola_set <- data.frame("label" = fola_set_df[,c("label")])
# convert it to a list 
fola_set_2 <- paste(fola_set$label, sep = ",")

# subset race data
fola_set_df_race <- data.frame("race" = metadata[,c("race")] )
rownames(fola_set_df_race) <- metadata$label

# ---- subset fola heatmap data ---- #
# reduce the white space in the heatmap but filtering columns where there is no data for foc
# first we extract only the foc rows using the same approach as for the metadata, but instead we perform on the binary matrix
fola_heat_df <- subset(binary_matrix, grepl("Fo._fsp._lactucae|Fo._fsp._matthiolae|Fo._Fo47", rownames(binary_matrix)))
# now we need to drop the empty columns 
fola_heat_df <- fola_heat_df[, colSums(fola_heat_df != 0) > 0]

# ---- Cluster the heatmap data (again)---- #
# normalisiation is mandatory for clustering, but as my data is binary - i did not normalise. 
# Compute hierarchical clustering of columns
fola_heat_df <- cluster_matrix(fola_heat_df, method ="ward.D2")

# ---- build the basic foc tree ---- #
#dop the non foc tip from the tree
tree_reduced <- drop.tip(tree, fola_set_2)
#build a tree from this data and add in out metadata
fola_tree <- ggtree(tree_reduced, ladderize = T ) %<+% metadata 

# ---- build full foc specific tree ---- #
# now make the tree pretty 
fola_tree_2 <- fola_tree +  
  geom_tiplab(aes(label = fsp), offset = 0.0003) +
  geom_tiplab(aes(label = isolate_code), color = "grey20", offset = 0.0018, linetype = "blank", geom = "text", align = TRUE) +
  geom_tippoint(aes(shape = source), size = 3) +
  geom_nodelab(geom='label', aes(label=label, subset= !is.na(as.numeric(label)) & as.numeric(label)> 40), nudge_x = -0.00039) +
  geom_rootedge() +
  theme(legend.position = "bottom")

# ---- add race data ---- #
#add extra scale so we can plot race with colour
fola_tree_3 <- fola_tree_2 + new_scale_fill()

fola_tree_4 <- gheatmap(fola_tree_3, fola_set_df_race, 
               offset = 0.0028, 
               width = 0.1,
               color = "grey",
               colnames = FALSE) +
  scale_fill_manual(name = "Race",
                    values = c("bisque", "gold4", "navy"), na.value = "white") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))

# ---- add candidate effector clusters ---- #

#add extra scale so we can plot race with colour
fola_tree_5 <- fola_tree_4 + new_scale_fill()

# add race data
fola_tree_6 <- gheatmap(fola_tree_5, fola_heat_df, offset=0.0035, colnames=T, colnames_angle=90, hjust=0.2, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 10)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "grey90", high = "black",
                        breaks = c("Absent","Present"),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())  +
   geom_treescale(x=0, y=.1, width = 0.004) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
    theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 8,
                                   colour = "white" )) +
  coord_cartesian(clip = "off") 

# plot the tree
plot(fola_tree_6)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/lactucae%20effector%20distib-1.png)<!-- -->

``` r
# save the tree
ggsave("HeatmapAndPhylo_LactucaeOnly.png", width = 20, height = 10)
```

### Misc

``` r
# ---- load session data --- #
session_data <- sessionInfo()
session_data
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] circlize_0.4.15       pheatmap_1.0.12       ggnewscale_0.4.9     
    ##  [4] RColorBrewer_1.1-3    textshape_1.7.3       ComplexHeatmap_2.15.4
    ##  [7] treeio_1.26.0         ggtreeExtra_1.13.0    ggtree_3.10.0        
    ## [10] phytools_2.1-1        maps_3.4.2            ape_5.7-1            
    ## [13] nortest_1.0-4         ggpubr_0.6.0          viridis_0.6.5        
    ## [16] viridisLite_0.4.2     ggthemes_5.0.0        lubridate_1.9.3      
    ## [19] forcats_1.0.0         stringr_1.5.1         purrr_1.0.2          
    ## [22] readr_2.1.5           tibble_3.2.1          ggplot2_3.4.4        
    ## [25] tidyverse_2.0.0       tidytree_0.4.6        tidyr_1.3.1          
    ## [28] dplyr_1.1.4          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] mnormt_2.1.1            gridExtra_2.3           phangorn_2.11.1        
    ##  [4] rlang_1.1.3             magrittr_2.0.3          clue_0.3-65            
    ##  [7] GetoptLong_1.0.5        matrixStats_1.2.0       compiler_4.3.1         
    ## [10] mgcv_1.9-1              systemfonts_1.0.5       png_0.1-8              
    ## [13] vctrs_0.6.5             combinat_0.0-8          quadprog_1.5-8         
    ## [16] shape_1.4.6             pkgconfig_2.0.3         crayon_1.5.2           
    ## [19] fastmap_1.1.1           magick_2.8.2            backports_1.4.1        
    ## [22] labeling_0.4.3          utf8_1.2.4              rmarkdown_2.25         
    ## [25] tzdb_0.4.0              ragg_1.2.7              xfun_0.41              
    ## [28] cachem_1.0.8            aplot_0.2.2             clusterGeneration_1.3.8
    ## [31] jsonlite_1.8.8          highr_0.10              cluster_2.1.6          
    ## [34] broom_1.0.5             parallel_4.3.1          R6_2.5.1               
    ## [37] stringi_1.8.3           car_3.1-2               numDeriv_2016.8-1.1    
    ## [40] Rcpp_1.0.12             iterators_1.0.14        knitr_1.45             
    ## [43] optimParallel_1.0-2     IRanges_2.36.0          splines_4.3.1          
    ## [46] Matrix_1.6-5            igraph_1.5.1            timechange_0.3.0       
    ## [49] tidyselect_1.2.0        rstudioapi_0.15.0       abind_1.4-5            
    ## [52] yaml_2.3.8              doParallel_1.0.17       codetools_0.2-19       
    ## [55] lattice_0.22-5          withr_3.0.0             coda_0.19-4            
    ## [58] evaluate_0.23           gridGraphics_0.5-1      pillar_1.9.0           
    ## [61] carData_3.0-5           stats4_4.3.1            foreach_1.5.2          
    ## [64] ggfun_0.1.4             generics_0.1.3          hms_1.1.3              
    ## [67] S4Vectors_0.40.2        munsell_0.5.0           scales_1.3.0           
    ## [70] glue_1.7.0              scatterplot3d_0.3-44    lazyeval_0.2.2         
    ## [73] tools_4.3.1             data.table_1.15.0       ggsignif_0.6.4         
    ## [76] fs_1.6.3                cowplot_1.1.3           fastmatch_1.1-4        
    ## [79] colorspace_2.1-0        nlme_3.1-164            patchwork_1.2.0        
    ## [82] cli_3.6.2               textshaping_0.3.7       fansi_1.0.6            
    ## [85] expm_0.999-9            gtable_0.3.4            rstatix_0.7.2          
    ## [88] yulab.utils_0.1.4       digest_0.6.34           BiocGenerics_0.48.1    
    ## [91] ggplotify_0.1.2         farver_2.1.1            rjson_0.2.21           
    ## [94] memoise_2.0.1           htmltools_0.5.7         lifecycle_1.0.4        
    ## [97] GlobalOptions_0.1.2     MASS_7.3-60.0.1
