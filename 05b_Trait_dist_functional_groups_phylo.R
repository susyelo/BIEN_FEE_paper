# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)
library(sp)
library(foreach)
library(raster)
library(viridis)
library(ggridges)
library(wesanderson)
library(lattice)


# Functions ---------------------------------------------------------------
source("./functions/check_functions.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/trait_distribution_functions.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

#2. Presence matrix of species
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv")
names(spPresence) = c("Species","Y","X")

#3. Richness raster as reference
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

#4. Biomes shapefiles
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

#5. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)

# Check raster and shapefiles are in the same projection
#proj4string(r_Total_Rich)==proj4string(biome_shp)

# Presence absence matrix of species in each biome -------------------------------------------------
TraitSpecies <- unique(Traits_phylo$species)
spMatrix_sub <- splistToMatrix(spPresence,TraitSpecies)

Biomes_pabs<-spMatrix_biomes(spMatrix = spMatrix_sub,
                             biome_shp = biome_shp,
                             raster_ref = r_Total_Rich)

# Compute distance matrix of trait between each pair of species  ----------
rownames(Traits_phylo)<-Traits_phylo$species

# Select traits to calculate the distances among species
traits<-c("Wood_density","Leaf_N","SLA","Seed_mass","Height","Leaf_P")


Dist_matrix<-compute_dist_matrix(Traits_phylo[,traits],metric="euclidean",
                                 center = TRUE,
                                 scale = TRUE) ## This can take a while



# Compute functional distinctiveness per biome ----------------------------

# Calculating relative abundance
Biome_relAbun<-make_relative(Biomes_pabs$Biomes_Abun)

Biomes_di = distinctiveness(Biome_relAbun, Dist_matrix)

Biomes_di_clean<-
  Biomes_di%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_di))%>%
  gather(key="Biome",value="Di",-species) %>% 
  filter(!is.na(Di))


biome_names=biome_shp$biomes

# Compute functional restrictiness per biome ----------------------------
## A 0 value indicates that the focal species is present in all the sites.
rest_species<-
  foreach(i=1:length(biome_names), .combine = rbind)%do%
  {
    indx<-which(rownames(Biomes_pabs$Biomes_pabs_cells)==biome_names[i])
    biome_PAbs_tmp<-Biomes_pabs$Biomes_pabs_cells[indx,]
    rest_species<-restrictedness(pres_matrix = biome_PAbs_tmp)
    rest_species$Biome<-biome_names[i]
    
    rest_species
  }

## Exclude all that are not present in any biome cell (which are not in Biome_Di)
# Usually species restricted to just one cell which is not classified in any of the biomes
rest_species<-
  rest_species %>% 
  filter(species%in%unique(Biomes_di_clean$species))


## Scaling Di values per biome
Biomes_di_clean<-
  Biomes_di_clean %>% 
  group_by(Biome) %>% 
  mutate(DiScale=rescaleRas01(Di))

## Merge Di and Ri
Biome_Di_Ri<-merge(Biomes_di_clean, rest_species)

Biome_Di_Ri$FunDi<-(Biome_Di_Ri$DiScale+Biome_Di_Ri$Ri)/2
Biome_Di_Ri$FunDi2<-Biome_Di_Ri$Di*Biome_Di_Ri$Ri

write.csv(Biome_Di_Ri, "./outputs/Biome_Di_Ri_phylo.csv")


## Heatmaps
# Total headmap


foreach (index=1:length(biome_names))%do%{
  
  png(paste("./figs/Di_Ri_heatmaps/Heatmap_", biome_names[index],".png",sep=""))
  print(Di_Ri_heatmaps(Biome_Di_Ri = Biome_Di_Ri, 
                 xvar = Biome_Di_Ri$Ri,
                 yvar = Biome_Di_Ri$DiScale,
                 Biome_toPlot = biome_names[index]))
  dev.off()
  
  
} 
