# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(foreach)

# functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. Presence of species in cells
## I am still not sure the order of the variables in the table
# TODO: Check this! 
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv",header = FALSE,
                     col.names = c("Species","Y","X"))


#2. Total_richness raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

## Include NA to the cells that have not information on them
r_Total_Rich[r_Total_Rich==0]<-NA

#spplot(r_Total_Rich)

# 3. Shapefiles
load("./data/base/Danilo_data/biomes_shp4susy.RData")

## Include cell number using the row and col numbers as reference
spPresence$cells<-cellFromRowCol(r_Total_Rich,spPresence$Y, spPresence$X)


# Extract species for each biome --------------------------------------------
# Ignoring montane TODO: do it with better polygons
biome_poly<-list(coniferous,dry,grasslands,
                 moist,prairies,savanna,taiga,
                 temperate.mixed,tropical.mixed,tundra)

names(biome_poly)<-c("coniferous","dry","grasslands",
                     "moist","prairies","savanna","taiga",
                     "temperate.mixed","tropical.mixed","tundra")

biome_richness<-foreach(i=1:length(biome_poly))  %do% {
  
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,biome_poly[[i]]))
  sp_list_tmp<-unique(spPresence$Species[spPresence$cells%in%cells_tmp])
  
}
names(biome_richness)<-names(biome_poly)


## Create similarity matrix
## Create a loop to calculate the similarity (number of species shared among biomes)

spSimilarity<-foreach(i=1:length(biome_richness), .combine='cbind') %:%
  foreach(j=1:length(biome_richness), .combine='c') %do% {
    length(intersect(biome_richness[[i]],biome_richness[[j]]))
  }

colnames(spSimilarity)<-names(biome_richness)
rownames(spSimilarity)<-names(biome_richness)

## Double check that the numbers are correct
#diag(spSimilarity)==unlist(lapply(biome_richness, n_distinct))


# Chordplot of similarities -----------------------------------------------


