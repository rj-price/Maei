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
tree_file <- "./MaeiTEFPhylo.treefile"
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

# build df for fsp so we can add it as a colour scale to the tree
fsp_df <- data.frame("fsp" = metadata[,c("fsp")] )
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
unrootedtree <- read.tree(tree_file)
# root the tree
tree <- root(unrootedtree, outgroup = c("F._graminearum_PH-1"))
# Build tree skeleton
p <- ggtree(tree,ladderize = F)  %<+% metadata

# ---- View the tree ---- #
# Useful for visualusing nodes etc 
p_nodes <- p + 
  geom_text2(aes(label = parent), hjust = -0.1, size = 3)+ # add node names
  geom_tiplab(aes(label = label), offset = 0.005) +
  coord_cartesian(clip = "off") # stop names being trimmed off
```

Now I have my basic tree, I can start to build something that will stand
alone. First I added the metadata (full name, the isolate code, and
race).

``` r
# ---- Build the tree plot ---- #

p2 <- p +  
  geom_treescale(x = 0, y = 1, width = 0.004) + 
  geom_tiplab(aes(label = full_name)) +
  geom_tiplab(aes(label = isolate_code), color = "black", offset = 0.015, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), color = "black", offset = 0.024, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source)) +
  geom_rootedge() +
  theme(legend.position = "bottom")

#add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()
# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.008, 
               width = 0.03,
               color = "black",
               colnames = FALSE) +
  scale_fill_manual(name = "Fsp",
                    values = c("blue","purple","goldenrod4","grey90","gold","brown", "lightpink","darkolivegreen3", "black", "tomato", "lavender", "tan", "palegreen4", "coral", "yellow"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

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
ggsave("BasicTEFPhylo.png", width = 30, height = 15)
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
  rename(isolate=isolate_code,
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
| F. oxysporum   | coriandrii   | 03-Feb       |   478 |                 315 |          65.4 |
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
  geom_histogram(fill = "#0c4c8a", colour = "black" ) +
  theme_bw()

# visualise the mimp distribution
cands_histo <- ggplot(stats_data, aes(x = candidate_effectors)) +
  geom_histogram(fill = "#0c4c8a", colour = "black" ) +
  theme_bw()

# visualise the mimp distribution
size_histo <- ggplot(stats_data, aes(x = assembly_size)) +
  geom_histogram(fill = "#0c4c8a", colour = "black" ) +
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
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw()

#Correlation between  the total number of effectors and Assembly Size?
#Visualize the relationship
effectors_v_assembly_size_relat <- ggplot(stats_data) +
  aes(x = candidate_effectors, y = assembly_size) +
  geom_point(colour = "#0c4c8a") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw()

#Correlation between the total number of mimps and Assembly Size?
#Visualize the relationship
mimps_v_assembly_size_relat <- ggplot(stats_data) +
  aes(x = mimps, y = assembly_size) +
  geom_point(colour = "#0c4c8a") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw()
```

``` r
# combine the plots 
sats_plots <- ggarrange(mimps_histo, cands_histo, size_histo, mimps_v_assembly_size_relat, effectors_v_assembly_size_relat, effectors_v_mimps_relat,
          ncol = 3, 
          nrow = 2)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
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
  rename(isolate=isolate_code,                                    
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
           colour="black",
           position= 'dodge',                          #Ensure the bars are not stacked. 
           stat='identity')+                           #Add the mimp or predicted effector content. 
  scale_fill_manual("Legend", values=c("candidate_effectors" = "darkolivegreen", "mimps" = "#DDE0DA"), label=c("mimps", "candidate effectors"))+ #for some reason the labels have to be written the other way round...?
  facet_grid(~species_group,
             scales = "free_x",                        # Let the x axis vary across facets.
             space = "free_x",                         # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+                            # Move the facet labels to the bottom.
  geom_point(aes(y=assembly_size /scale_right,         #Plot assemble size over the top of the bar chart.
                 colour = "Assembly Size",             #Add assembly size to the legend. 
                 group = 1))+ 
  scale_colour_manual(" ", values=c("Assembly Size" = "black"))+
  theme_bw()+
  theme(legend.box="verticle",
        legend.title = element_blank())+
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90,         #Adjust the text orientation on the x axis
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 8))+
  theme(axis.text = element_text(size = 8))+
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
ggsave("StatsOverview.png", width = 15, height = 10)
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

I plotted the mimp and effector content from the metadata and visualised
it using `ggplot.`

``` r
# ---- Prepare Data for plotting --- #
#Extract the isolate, assembly size, total number of mimps and effectors coloumns.
race_plot_data <- select(metadata,"species", "species_group", "fsp", "race" ,"isolate_code","genome_size","no._mimps","no._cand_effs") %>%
#Rename the columns to reduce the long titles
  rename(isolate=isolate_code,                                    
         mimps=no._mimps,
         candidate_effectors=no._cand_effs,
         assembly_size =genome_size) %>%  
#We need to drop rows which were not included in the Maei analysis
  drop_na(candidate_effectors) %>%
#Merge/group the mimps and effector columns so that both can be plotted per strain/isolate.
  pivot_longer(cols = c(mimps,candidate_effectors), names_to="Legend", values_to="mimps_and_candidate_effectors") %>%  
  mutate(Legend = factor(Legend, levels=c('mimps','candidate_effectors')))


# ---- subset metadata ---- #
# extract only the fsp we are interested in 
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
       y = "Count") +
  scale_fill_manual(values = c("candidate_effectors" = "darkolivegreen", "mimps" = "#DDE0DA"), label=c("mimps", "candidate effectors")) +
  theme(strip.text.y = element_blank(),  # remove the side names as we have this shown in colour now. 
        panel.grid.major = element_line(), # put the lines back in 
        legend.position = "bottom", 
        legend.title = element_blank())
#plot it 
plot(mimpsandcandeffs)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/effector%20and%20mimp%20race%20plots-1.png)<!-- -->

``` r
#save the plot
ggsave("MimpsAndCandEffs_CubenseOnly.png", width = 20, height = 10)
```

## Candidate effector distribution

Now I have an understanding of the effector distribution and the
phylogeny, I combined the phylogenies and effector profiles to generate
a heat map. For this, I can use the tree (p4) already generated and add
the heatmap data I loaded initially (data).

First, I need to prepare the heatmap data. As we are looking at
clustered (0.65% ID, cd-hit (v…)) and filtered sequences (SignalP
(v5.06) and EffectorP (v2.0.1) extracted from BLAST hits, the total
number of candidate effectors per cluster per isolate varies, but I want
to just look at Presence/. In order to do this, I converted the heatmap
data matrix to a binary data frame.

``` r
# ---- Prep the heatmap data ---- #
# remove the row names temporarily 
rownames_mat<- data[,1]
mat_data<- as.matrix(data[,-1])
#make data frame binary 
binary_matrix <- as.matrix(mat_data)
binary_matrix[binary_matrix > 0] <- 1
#put rownames back
rownames(binary_matrix)<-rownames_mat
```

Next, I clustered the data in the binary data matrix, so that it will be
ordered when I visualise the candidate effector heatmap.

``` r
# ---- Cluster the heatmap data ---- #
# normalisiation is mandatory for clustering, but as my data is binary - i did not normalise. 
# Compute hierarchical clustering of columns
heatmap_dat <- cluster_matrix(binary_matrix, dim = 'col', method ="ward.D2")
```

Next, in theory, I add the heatmap to the tree already built… but I cant
get the offsets to align consistently and it overlaps if i just add
heatmap_dat to p4, so I have to rebuild p4.

``` r
# ---- Build the full figure ---- #
p2 <- p +  
  geom_treescale(x = 0, y = 1, width = 0.004) + 
  geom_tiplab(aes(label = full_name), offset = 0.001) +
  geom_tiplab(aes(label = isolate_code), color = "black", offset = 0.028, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), color = "black", offset = 0.047, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source), size = 3) +
  geom_rootedge() +
  theme(legend.position = "bottom")

#add extra scale so we can plot fsp with colour
p3 <- p2 + new_scale_fill()
# add race data
p4 <- gheatmap(p3, fsp_df,
               offset = 0.02, 
               width = 0.03,
               color = "black",
               colnames = FALSE) +
  scale_fill_manual(name = "Fsp",
                    values = c("blue","purple","goldenrod4","grey90","gold","brown", "lightpink","darkolivegreen3", "black", "tomato", "lavender", "tan", "palegreen4", "coral", "yellow"), na.value = "grey") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "vertical", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
# add extra discrete scale
p5 <- p4 + new_scale_fill()

# add effector heatmap
p6 <-gheatmap(p5, heatmap_dat, offset=0.07, colnames=FALSE, legend_title="Presence/\nAbsence", color = NULL,  width = 1.5)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "white", high = "darkolivegreen",
                        breaks = c("Absent","Present"),
                        na.value = "grey")+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box = "horizontal", legend.margin = margin())+
  guides(shape = guide_legend(override.aes = list(size = 2)))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
plot(p6)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/final%20heatmap-1.png)<!-- -->

``` r
# save the output
ggsave("HeatmapAndPhylo.png", width = 30, height = 10)
```

### Candidate effector distribution in Fo. fsp. cubense

Now I have my overall heatmap, I want to look at some of the fsp in more
detail - particularly those for which we have multiple races available.
I subset my metadata, searched for all caseses that did not match the
regex “Fo.\_fsp.\_cubense” and dropped them from the tree using the
`drop.tip` function. I the reconstructed my tree using just the Foc tef
phylo and foc effector profiles.

``` r
# ---- subset foc metadata ---- #
#identify all rows in the metadata which do not contain cubense 
foc_set_df <- subset(metadata, !grepl("Fo._fsp._cubense|Fo._Fo47", label))
#subset just tip labels
foc_set <- data.frame("label" = foc_set_df[,c("label")])
# convert it to a list 
foc_set_2 <- paste(foc_set$label, sep = ",")

# ---- subset foc heatmap data ---- #
# reduce the white space in the heatmap but filtering columns where there is no data for foc
# first we extract only the foc rows using the same approach as for the metadata, but instead we perform on the binary matrix
foc_heat_df <- subset(binary_matrix, grepl("Fo._fsp._cubense|Fo._Fo47", rownames(binary_matrix)))
# now we need to drop the empty columns 
foc_heat_df <- foc_heat_df[, colSums(foc_heat_df != 0) > 0]

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
  #geom_tiplab(aes(label = full_name), offset = 0.00001) +
  geom_tiplab(aes(label = isolate_code), color = "black", offset = 0.0003, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), offset = 0.0038, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source), size = 3) +
  geom_rootedge() +
  theme(legend.position = "bottom")

#add extra scale so we can plot fsp with colour
foc_tree_3 <- foc_tree_2 + new_scale_fill()
# add race data
foc_tree_4 <- gheatmap(foc_tree_3, foc_heat_df, offset=0.008, colnames=T, colnames_angle=90, hjust=1, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 4)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "white", high = "darkolivegreen",
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
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
# plot the tree
plot(foc_tree_4)
```

![](AnalysisCandidateEffectorSets_files/figure-gfm/cubense%20effector%20distib-1.png)<!-- -->

``` r
# save the tree
ggsave("HeatmapAndPhylo_CubenseOnly.png", width = 20, height = 10)
```

We can see that the effector profiles differ based on race and overall
tef1 phylogeny. Interestingly there is one candidate which is shared
across the R1 and suspected STR4 isolates but not found in the TR4
isolates - cognate R gene?

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

### Candidate effector distribution in Fo. fsp. lactucae and Fo. fsp. apii

``` r
# ---- subset metadata ---- #
#identify all rows in the metadata which do not contain apii and coriandrii 
foa_c_set_df <- subset(metadata, !grepl("Fo._fsp._coriandrii|Fo._fsp._apii|Fo._Fo47", label))
#subset just tip labels
foa_c_set <- data.frame("label" = foa_c_set_df[,c("label")])
# convert it to a list 
foa_c_set_2 <- paste(foa_c_set$label, sep = ",")

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
  geom_tiplab(aes(label = fsp), color = "black", offset = 0.0003, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = isolate_code), color = "black", offset = 0.002, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), offset = 0.0042, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source), size = 3) +
  geom_rootedge() +
  theme(legend.position = "bottom")

#add extra scale so we can plot fsp with colour
foa_c_tree_3 <- foa_c_tree_2 + new_scale_fill()
# add race data
foa_c_tree_4 <- gheatmap(foa_c_tree_3, foa_c_heat_df, offset=0.0065, colnames=T, colnames_angle=90, hjust=1, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 4)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "white", high = "darkolivegreen",
                        breaks = c("Absent","Present"),
                        na.value = "grey")+
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
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
# plot the tree
plot(foa_c_tree_4)
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
  geom_tiplab(aes(label = isolate_code), color = "black", offset = 0.002, linetype = "blank", geom = "text", align = TRUE) +
  geom_tiplab(aes(label = race), offset = 0.0042, linetype = "blank", geom = "text", align = TRUE)+
  geom_tippoint(aes(shape = source), size = 3) +
  geom_rootedge() +
  theme(legend.position = "bottom")

#add extra scale so we can plot fsp with colour
fola_tree_3 <- fola_tree_2 + new_scale_fill()

# add race data
fola_tree_4 <- gheatmap(fola_tree_3, fola_heat_df, offset=0.0065, colnames=T, colnames_angle=90, hjust=1, font.size=3, legend_title="Presence/\nAbsence", color = "grey",  width = 4)  +
  scale_fill_continuous(name = "Presence/\nAbsence",
                        low = "white", high = "darkolivegreen",
                        breaks = c("Absent","Present"),
                        na.value = "grey") +
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
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
# plot the tree
plot(fola_tree_4)
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.12    ggnewscale_0.4.9   RColorBrewer_1.1-3 textshape_1.7.3   
    ##  [5] ggtreeExtra_1.13.0 ggtree_3.10.0      phytools_2.1-1     maps_3.4.2        
    ##  [9] ape_5.7-1          nortest_1.0-4      ggpubr_0.6.0       viridis_0.6.5     
    ## [13] viridisLite_0.4.2  ggthemes_5.0.0     lubridate_1.9.3    forcats_1.0.0     
    ## [17] stringr_1.5.1      purrr_1.0.2        readr_2.1.5        tibble_3.2.1      
    ## [21] ggplot2_3.4.4      tidyverse_2.0.0    tidytree_0.4.6     tidyr_1.3.1       
    ## [25] dplyr_1.1.4       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] mnormt_2.1.1            gridExtra_2.3           phangorn_2.11.1        
    ##  [4] rlang_1.1.3             magrittr_2.0.3          compiler_4.3.1         
    ##  [7] mgcv_1.9-1              systemfonts_1.0.5       vctrs_0.6.5            
    ## [10] combinat_0.0-8          quadprog_1.5-8          pkgconfig_2.0.3        
    ## [13] fastmap_1.1.1           backports_1.4.1         labeling_0.4.3         
    ## [16] utf8_1.2.4              rmarkdown_2.25          tzdb_0.4.0             
    ## [19] ragg_1.2.7              xfun_0.41               cachem_1.0.8           
    ## [22] aplot_0.2.2             clusterGeneration_1.3.8 jsonlite_1.8.8         
    ## [25] highr_0.10              broom_1.0.5             parallel_4.3.1         
    ## [28] R6_2.5.1                stringi_1.8.3           car_3.1-2              
    ## [31] numDeriv_2016.8-1.1     Rcpp_1.0.12             iterators_1.0.14       
    ## [34] knitr_1.45              optimParallel_1.0-2     splines_4.3.1          
    ## [37] Matrix_1.6-5            igraph_1.5.1            timechange_0.3.0       
    ## [40] tidyselect_1.2.0        rstudioapi_0.15.0       abind_1.4-5            
    ## [43] yaml_2.3.8              doParallel_1.0.17       codetools_0.2-19       
    ## [46] lattice_0.22-5          treeio_1.26.0           withr_3.0.0            
    ## [49] coda_0.19-4             evaluate_0.23           gridGraphics_0.5-1     
    ## [52] pillar_1.9.0            carData_3.0-5           foreach_1.5.2          
    ## [55] ggfun_0.1.4             generics_0.1.3          hms_1.1.3              
    ## [58] munsell_0.5.0           scales_1.3.0            glue_1.7.0             
    ## [61] scatterplot3d_0.3-44    lazyeval_0.2.2          tools_4.3.1            
    ## [64] data.table_1.15.0       ggsignif_0.6.4          fs_1.6.3               
    ## [67] cowplot_1.1.3           fastmatch_1.1-4         grid_4.3.1             
    ## [70] colorspace_2.1-0        nlme_3.1-164            patchwork_1.2.0        
    ## [73] cli_3.6.2               textshaping_0.3.7       fansi_1.0.6            
    ## [76] expm_0.999-9            gtable_0.3.4            rstatix_0.7.2          
    ## [79] yulab.utils_0.1.4       digest_0.6.34           ggplotify_0.1.2        
    ## [82] farver_2.1.1            memoise_2.0.1           htmltools_0.5.7        
    ## [85] lifecycle_1.0.4         MASS_7.3-60.0.1
