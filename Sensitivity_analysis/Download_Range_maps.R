# Libraries ---------------------------------------------------------------
library(BIEN)
library(ape)
library(picante)
library(tidyverse)
library(raster)


# Data --------------------------------------------------------------------
# 1. Traits data
Trait_BIEN_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv", row.names=1)
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ","_",Trait_BIEN_df$scrubbed_species_binomial)

#2. Phylogenetic data
Seed_phylo<-read.tree("./data/base/big_seed_plant_trees_v0.1/ALLMB.tre")


#3. Baseline raster
r<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")


# Subset data -------------------------------------------------------------
Trait_BIEN_df <-
Trait_BIEN_df %>% 
  filter(scrubbed_species_binomial%in%Seed_phylo$tip.label)


# Download range maps -----------------------------------------------------
range_maps<-BIEN_ranges_load_species(Trait_BIEN_df$scrubbed_species_binomial)
saveRDS(range_maps, "./data/processed/BIEN3.0SeedPhylo_maps.rds")

range_maps_df<-spTransform(range_maps, CRS("+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
+ellps=WGS84 +towgs84=0,0,0")) 

resol<-2
r<-raster(xmn = -167, xmx = -29,ymn = -57, ymx = 83,
          crs = CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
res(r) <- resol

p <- rasterToPolygons(r)
p$cell <- 1:ncell(r)
p$layer <- NULL

v <- intersect(range_maps, p)

m <- as.matrix(table(v$cell, v$species))
m[1:10, ]


range_maps_ras<- rasterize(range_maps, r)

