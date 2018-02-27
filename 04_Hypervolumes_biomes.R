# libraries ---------------------------------------------------------------
library(hypervolumes)
library(raster)
library(rworldmap)

# data --------------------------------------------------------------------
# Olson shapefiles and raster files
load("./data/base/Danilo_data/biomes_shp4susy.RData")
load("./data/base/Olson_raster.RData")

## Presence/absence matrix of species with traits
load("./data/processed/PresAbs_matrix_TraitRange.RData")

#Richness raster of the species with traits
richness_ras<-raster("./data/processed/TraitRichness_raster.tif")
