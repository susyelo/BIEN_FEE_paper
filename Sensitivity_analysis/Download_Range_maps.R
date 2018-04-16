# Libraries ---------------------------------------------------------------
library(BIEN)
library(ape)
library(picante)
library(tidyverse)
library(raster)
library(foreach)


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
saveRDS(range_maps, "./data/base/BIEN3.0SeedPhylo_maps.rds")



## Create baseline raster: to extract the presence/absence of species in each cell
# Base for resolution and projection
resol<-2
r_base<-raster(xmn = -167, xmx = -29,ymn = -57, ymx = 83,
               crs = CRS(proj4string(r)))

res(r_base) <- resol
values(r_base)<-1:ncell(r_base)
r_base<-projectRaster(r_base, r)

p <- rasterToPolygons(r_base)
p$cell <- 1:ncell(r_base)
p$layer <- NULL

# Absence/presence matrix of species
num<-seq(1,length(range_maps),100)
tmp_num<-data.frame(n1=num, n2=c(num[-1]+1,length(range_maps)))

Cells_sp<-
foreach(i=1:nrow(tmp_num), .combine = rbind)%do%{
  
  col1<-tmp_num$n1[i]
  col2<-tmp_num$n2[i]
  print(paste("Extract ranges from",range_maps$species[col1], "to", range_maps$species[col2]))
  range_maps_sp<-range_maps[c(col1:col2),]
  range_maps_df<-spTransform(range_maps_sp, CRS(proj4string(r)))
  
  v <- over(range_maps_df,p)
  v$species<-range_maps_df$species
  
  v
}

