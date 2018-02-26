# Libraries ---------------------------------------------------------------
library(letsR)
library(tidyverse)
library(gridExtra)

# Data --------------------------------------------------------------------
# trait data
trait_df<-read.csv("./outputs/BIEN_trait_GrowthForm.csv")

# Occurrence data
occ_df<-readRDS("./data/base/2018_02_08_BIEN_traitOccurData.RData")


# Delete names spaces
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ", "_",Trait_BIEN_df$scrubbed_species_binomial)

# Filter Occurrence data --------------------------------------------------
Occ_BIEN_sub<-
  occ_df %>% 
  filter(scrubbed_species_binomial %in% unique(trait_df$scrubbed_species_binomial))


# Create species richness map
xy<-data.frame(x=Occ_BIEN_sub$longitude, y=Occ_BIEN_sub$latitude)

# Presence/abscence raster of Solanum ---
BIEN_trait_grids<-lets.presab.points(xy, Occ_BIEN_sub$scrubbed_species_binomial, 
                                     resol=1,count = TRUE,
                                     xmn = -167, xmx = -29,ymn = -57, ymx = 83)

save(BIEN_trait_grids, file="./outputs/BIEN_traits_PresAbs1degree.RData")

BIEN_trait_grids2d<-lets.presab.points(xy, Occ_BIEN_sub$scrubbed_species_binomial, 
                                     resol=2,count = TRUE,
                                     xmn = -167, xmx = -29,ymn = -57, ymx = 83)


BIEN_trait_grids2d$Richness_Raster[BIEN_trait_grids2d$Richness_Raster==0]<-NA

plot(log10(BIEN_trait_grids$Richness_Raster))


# Trait maps --------------------------------------------------------------
# range size
range_size<-lets.rangesize(BIEN_trait_grids, units="squaremeter")

Range_size_map<- lets.maplizer(BIEN_trait_grids,
                               range_size,
                               BIEN_trait_grids$Species_name,
                               ras = TRUE)

# Height
trait_df_wide<-
  trait_df  %>% 
  spread(trait_name,trait_value) %>% 
  filter(scrubbed_species_binomial%in%BIEN_trait_grids$Species_name)

Height<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(height=mean(whole_plant_height,na.rm=TRUE))


Height_map<- lets.maplizer(BIEN_trait_grids,
                           log(Height$height),
                           BIEN_trait_grids$Species_name,
                           func = mean,
                           ras = TRUE)

plot(Height_map$Raster)


# Specific leaf area
## really struggling to get this trait work (no enough data to calculate it)

SLA<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(SLA=mean(whole_plant_leaf_area_per_whole_plant_leaf_dry_mass,na.rm=TRUE))

SLA_map<- lets.maplizer(BIEN_trait_grids,
                        SLA$SLA,
                        BIEN_trait_grids$Species_name,
                        func = mean,
                        ras = TRUE)

plot(SLA_map$Raster)


# Seed mass
Seed_mass<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(seed_mass=mean(seed_mass,na.rm=TRUE))

Seed_mass_map<- lets.maplizer(BIEN_trait_grids,
                              log(Seed_mass$seed_mass),
                              BIEN_trait_grids$Species_name,
                              func = mean,
                              ras = TRUE)
plot(Seed_mass_map$Raster,col=terrain.colors(100))

spplot(Seed_mass_map$Raster)


# Leaf N
Leaf_N<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(Leaf_N=mean(leaf_nitrogen_content_per_leaf_dry_mass,na.rm=TRUE))

Leaf_N_map<- lets.maplizer(BIEN_trait_grids,
                           Leaf_N$Leaf_N,
                           BIEN_trait_grids$Species_name,
                           func = mean,
                           ras = TRUE)

plot(Leaf_N_map$Raster,col=terrain.colors(100))
spplot(Leaf_N_map$Raster)

# Leaf P
Leaf_P<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(Leaf_P=mean(leaf_phosphorus_content_per_leaf_dry_mass,na.rm=TRUE))

Leaf_P_map<- lets.maplizer(BIEN_trait_grids,
                           Leaf_P$Leaf_P,
                           BIEN_trait_grids$Species_name,
                           func = mean,
                           ras = TRUE)

plot(Leaf_P_map$Raster,col=terrain.colors(100))
spplot(Leaf_P_map$Raster)

# Wood density
Wood_dens<-  
  trait_df_wide %>% 
  group_by(scrubbed_species_binomial) %>% 
  summarise(Wood_dens=mean(whole_plant_woodiness,na.rm=TRUE))

Wood_dens_map<- lets.maplizer(BIEN_trait_grids,
                           Wood_dens$Wood_dens,
                           BIEN_trait_grids$Species_name,
                           func = mean,
                           ras = TRUE)

Wood_dens_map_sd<- lets.maplizer(BIEN_trait_grids,
                              Wood_dens$Wood_dens,
                              BIEN_trait_grids$Species_name,
                              func = sd,
                              ras = TRUE)

plot(Wood_dens_map$Raster,col=terrain.colors(100))

spplot(Wood_dens_map_sd$Raster)

spplot(Wood_dens_map$Raster)
