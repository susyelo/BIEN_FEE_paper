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
### For loops for each species in a 200x200km resolution map
#r_base<-projectRaster(r, 
#                      crs=crs(r),
#                      res=2e+05)
r_base<-r
values(r_base)<-1:ncell(r_base)
names(r_base)<-"cell"

## Convert raster to cell polygons
p <- rasterToPolygons(r_base)

n_species<-length(range_maps$species)

Cells_sp<-
  foreach(i=1:n_species, .combine = rbind)%do%{
    
    print(paste("Extract ranges from",range_maps$species[i]))
    range_maps_sp<-range_maps[i,]
    range_maps_df<-spTransform(range_maps_sp, CRS(proj4string(r)))
    
    v <-as.data.frame(over(range_maps_df,p,returnList = TRUE))
    v$species<-range_maps_df$species
    
    v
  }

#saveRDS(Cells_sp, "outputs/Cells_sp_BIEN_200km.rds")
saveRDS(Cells_sp, "outputs/Cells_sp_BIEN_100km.rds")

pb_matrix<-as.matrix(table(Cells_sp$cell,Cells_sp$species))

cell_richness<-data.frame(N_sp=as.numeric(rowSums(pb_matrix)),
                          Cell=as.numeric(rownames(pb_matrix)))


## Species richness with all the species with traits that are in the plant seed phylogeny
values(r_base)<-NA
values(r_base)[cell_richness$Cell]<-cell_richness$N_sp
spplot(r_base)


## Total proportion of species with trait information in the study (i.e., species with traits and included in the Smith and Brown phylogeny)
## Change resolution and aggregate the richness from 100x100km to 200x200km. 


r200<-aggregate(r, fact=2, fun=sum, expand=FALSE)## This calculates the average values of number of species in 4 cells

r_base<-crop(r_base, extent(r200))

indx<-which(values(!is.na(r200)))

values(r_base)[indx][is.na(values(r_base)[indx])]<-0


## Proportion of species with trait info
Total_prop_map<-r_base/r200
spplot(Total_prop_map)
