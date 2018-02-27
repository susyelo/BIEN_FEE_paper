# libraries ---------------------------------------------------------------
library(hypervolumes)
library(raster)
library(rworldmap)
library(tidyverse)

# data --------------------------------------------------------------------
# Olson shapefiles and raster files
load("./data/base/Danilo_data/biomes_shp4susy.RData")


## Presence/absence matrix of species with traits
load("./data/processed/PresAbs_matrix_TraitRange.RData")

#Richness raster of the species with traits
richness_ras<-raster("./data/processed/TraitRichness_raster.tif")

biomes_mean_PCs<-read.table("./data/base/Danilo_data/biomes.csv")
load("./data/base/Danilo_data/biomes_shp4susy.RData")



# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

## Find cells in raster using coordinates as 
cells <- cellFromXY(richness_ras, biomes_mean_PCs[,c("Long10","Lat10")])

cells_dry<-cellFromPolygon(richness_ras,dry)

## To check if the cells correspond to the biome
r_tmp<-richness_ras
r_tmp[cells_dry[[1]]]<-500000
plot(r_tmp)


# Extract list of species for each biome
tmp_matrix<-as.data.frame(spMatrix[cells_dry[[1]],])

tmp_matrix %>% 
  gather(key="Species", value = "Presence") %>% 
  filter(Presence!=0) %>% 
  select(Species) %>% 
  unique

