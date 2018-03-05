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
                     col.names = c("Species","X","Y"))

#2. Total_richness raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

## Include NA to the cells that have not information on them
r_Total_Rich[r_Total_Rich==0]<-NA

#spplot(r_Total_Rich)

# 3. Shapefiles
load("./data/Danilo_data/biomes_shp4susy.RData")

# Extract species for each biome --------------------------------------------
# Ignoring montane TODO: do it with better polygons
biome_poly<-list(coniferous,dry,grasslands,
                 moist,prairies,savanna,taiga,
                 temperate.mixed,tropical.mixed,tundra)

names(biome_poly)<-c("coniferous","dry","grasslands",
                     "moist","prairies","savanna","taiga",
                     "temperate.mixed","tropical.mixed","tundra")

biome_richness<-foreach(i=1:length(biome_poly))  %do% {
  
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,coniferous))
  sp_list_tmp<-spFromCell(d=spPresence,cell = cells_tmp, r=r_Total_Rich)
  
}

names(biome_richness)<-names(biome_poly)


## Create similarity matrix

# Make matrices to store the final results, temporary sums and counts
count.matrix<-matrix(0,nrow=length(biome_richness),ncol=length(biome_richness))

## Create a loop to calculate the similarity
length(intersect(biome_richness$coniferous,biome_richness$dry))


  


