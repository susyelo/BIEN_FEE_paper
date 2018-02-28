# libraries ---------------------------------------------------------------
library(hypervolume)
library(raster)
library(rworldmap)
library(tidyverse)
library(visdat)
library(RColorBrewer)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. Olson shapefiles and raster files
load("./data/base/Danilo_data/biomes_shp4susy.RData")

## 2. Presence/absence matrix of species with traits
load("./data/processed/PresAbs_matrix_TraitRange.RData")

# 3. Richness raster of the species with traits
richness_ras<-raster("./data/processed/TraitRichness_raster.tif")

# 4. Biomes polygons
load("./data/base/Danilo_data/biomes_shp4susy.RData")

# 5. Trait data
trait_BIEN<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv")
trait_BIEN$scrubbed_species_binomial<-gsub(" ","_",trait_BIEN$scrubbed_species_binomial)



# Explore data from traits ------------------------------------------------
vis_dat(trait_BIEN, 
        sort_type = FALSE)

vis_miss(trait_BIEN)

vis_miss(trait_BIEN,
         sort_miss = TRUE)

# Extract species list per biome ------------------------------------------

## Find cells in raster using coordinates as 
cells_dry<-cellFromPolygon(richness_ras,coniferous)

## To check if the cells correspond to the biome
r_tmp<-richness_ras
r_tmp[cells_dry[[1]]]<-500000
plot(r_tmp)

# Extract list of species
biome_matrix<-as.data.frame(spMatrix[cells_dry[[1]],])

sp_list<-biome_matrix %>% 
  gather(key="Species", value = "Presence") %>% 
  filter(Presence!=0) %>% 
  dplyr::select(Species) %>% 
  unique


# Calculate hypervolume for the whole biome -------------------------------
var_names<-c("stem_wood_density",
             "whole_plant_leaf_area_per_whole_plant_leaf_dry_mass",
             "seed_mass",
             "whole_plant_height")

### Hypervolume for the species list
trait_sp_biomes<-trait_BIEN %>% 
  filter(scrubbed_species_binomial%in%sp_list$Species) %>% 
  select(var_names)


## I cannot use hypervolumes with NAs
# Without NA
tmp<-na.omit(trait_sp_biomes)
tmp<-scale(tmp[,1:ncol(tmp)])

hyper_tmp<-hypervolume(tmp,name="dry")

hyper_gaussian<-hypervolume_gaussian(tmp,name="dry")

hyper_gaussian_moist<-hypervolume_gaussian(tmp,name="moist")

hyper_gaussian_conifer<-hypervolume_gaussian(tmp,name="conifers")


##plotting hypevolumes
#Log transform amphibian hypervolume
plot(hyper_gaussian_moist)
plot(hyper_gaussian_moist,show.3d=TRUE,plot.3d.axes.id=2:4,cex.random=3,cex.data=6,
     show.legend=TRUE,point.alpha.min=0.5,point.dark.factor=1)

plot(hyper_gaussian,show.3d=TRUE,plot.3d.axes.id=2:4,cex.random=3,cex.data=6,
     show.legend=TRUE,point.alpha.min=0.5,point.dark.factor=1)

#Plotting all three amniote hypervolumes together
#Log transformed hypervolumes
plot(hypervolume_join(hyper_gaussian_moist,hyper_gaussian,hyper_gaussian_conifer),num.points.max.random=6000,contour.lwd=1.5,colors=c(brewer.pal(n=3,"Set1"))
     ,show.legend=FALSE)
legend("bottomleft",legend = c("moist","dry","conifer"),text.col=c(brewer.pal(n=3,"Set1")),bty="n",cex=1.1,text.font=2)
