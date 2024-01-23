#Jamie Pike 26/3/2021

#Build graph for mimp and effector distribution relative to genome, and test for correlation. 

#Set libraries
#-------------

#install.packages("ggplo2")
#install.packages("tidyr")
library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)
library(forcats)

#Load and Prepare Data 
#---------------------

#Set working directory 
setwd("/Volumes/Jamie_EXT/Fusarium_data/Maei_Analysis/Cubense_and_Species_Complex/Statistics")
#Read the csv file in as a table
StatsCSV <- read.table("MimpAndEffectorStats.csv", header = T, sep=",")

StatsData <- select(StatsCSV,"Strain","Assembly.Size..Mb.","Total.number.of.predicted.mimps","Total.number.of.predicted.effectors") %>%
  #Rename the columns to reduce the long titles
  rename(Isolate=Strain, 
         Mimps=Total.number.of.predicted.mimps,
         Effectors=Total.number.of.predicted.effectors,
         Assembly_Size = Assembly.Size..Mb.)

#Perform correlation statistics
#------------------------------
#Correlation coefficient calculated using the cor() function. The Pearson correlation is computed by default with the cor() function.

#Correlation between the total number of effectors and  the total number of mimps?
#Visualize the relationship
ggplot(StatsData) +
  aes(x = Effectors, y = Mimps) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()
#Correlation test:
EffectorsVsMimps <- cor.test(StatsData$Effectors, StatsData$Mimps) 
EffectorsVsMimps

#Correlation between  the total number of effectors and Assembly Size?
#Visualize the relationship
ggplot(StatsData) +
  aes(x = Effectors, y = Assembly_Size) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()
#Correlation test:
EffectorsVsAssemblySize <- cor.test(StatsData$Effectors, StatsData$Assembly_Size) 
EffectorsVsAssemblySize

#Correlation between the total number of mimps and Assembly Size?
#Visualize the relationship
ggplot(StatsData) +
  aes(x = Mimps, y = Assembly_Size) +
  geom_point(colour = "#0c4c8a") +
  theme_minimal()
#Correlation test:
MimpsVsAssemblySize <- cor.test(StatsData$Mimps, StatsData$Assembly_Size)
MimpsVsAssemblySize

#Note that the correlation between variables X and Y is equal to the correlation between variables Y and X so the order of the variables in the cor() function does not matter.



#Prepare Data for plotting 
#-------------------------

#Extract the strain, assembly size, total number of mimps and effectors coloumns.ß
StatsData <- select(StatsCSV,"Strain","Assembly.Size..Mb.","Total.number.of.predicted.mimps","Total.number.of.predicted.effectors") %>%
  #Rename the columns to reduce the long titles
  rename(Isolate=Strain, 
         Mimps=Total.number.of.predicted.mimps,
         Effectors=Total.number.of.predicted.effectors,
         Assembly_Size = Assembly.Size..Mb.) %>%
  #Merge/group the mimps and effector columns so that both can be plotted per strain/isolate.
  pivot_longer(cols = c(Mimps,Effectors), names_to="Gene_type", values_to="Mimps_and_Effectors") %>%
  mutate(Gene_type = factor(Gene_type, levels=c('Mimps','Effectors')))

##Check the correct columns are extracted.
print(StatsData)


#Build plot
#----------

#Generate scale for Assembly size data 
scaleRight <- 100 / max(StatsData$Mimps_and_Effectors)

#Build Plot
ggplot(aes(x=reorder(Isolate, Mimps_and_Effectors)), data = StatsData)+ #Create X axis, which contains all strains/isolates assessed. 
  geom_bar(aes(y=Mimps_and_Effectors, fill = Gene_type), position= 'dodge', stat='identity')+ #Add the mimp or predicted effector content. 
  geom_point(aes(y=Assembly_Size /scaleRight, colour = "Assembly Size", group = 1))+ #Plot assemblt size over the top of the bar chart.
  scale_colour_manual(" ", values=c("Assembly Size" = "black"))+
  # scale_fill_manual("",values="red")+
  theme(legend.box="verticle")+
  scale_fill_discrete(name = "Legend")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ #Adjust the text orentation on the x axis
  theme(panel.background = element_blank(), panel.border = element_blank())+ #Remove the grey background 
  xlab("Isolate")+
  scale_y_continuous(name= "Total number of predicted\nmimps and effector hits", breaks = scales::pretty_breaks(n = 20), 
                     sec.axis = sec_axis( trans=~./1*scaleRight, name="Size of Assembly (Mb)", breaks = scales::pretty_breaks(n = 20))) #Increase ticks on y axis.

                     