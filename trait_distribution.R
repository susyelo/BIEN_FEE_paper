# Libraries ---------------------------------------------------------------
library(tidyverse)
library(raster)
library(foreach)
library(plyr)
library(ggridges)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. List of species per biome (all species)
load("./outputs/Biome_list_species.RData")

# 2. Trait data BIEN with growth form
trait_BIEN<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv")
trait_BIEN$scrubbed_species_binomial<-gsub(" ","_",trait_BIEN$scrubbed_species_binomial)

# 3. Biomes shapefile
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")



# Filter species with traits ----------------------------------------------
indx<-which(biome_richness$Coniferous_Forests%in%trait_BIEN$scrubbed_species_binomial)
n_distinct(biome_richness$Coniferous_Forests[indx])

biome_richness_traits<-lapply(biome_richness, 
                              function(x)
                                x[which(x%in%trait_BIEN$scrubbed_species_binomial)])

#Convert values in dataframe 
biome_sp_trait_df <- ldply(biome_richness_traits, data.frame)
colnames(biome_sp_trait_df)<-c("Biome","Species")

# Merge trait and biome dataframes -----------------------------------------------
DF_trait_biome<-merge(trait_BIEN,biome_sp_trait_df,by.x="scrubbed_species_binomial",by.y="Species")

DF_trait_biome$Biome<-factor(DF_trait_biome$Biome,levels = c("Moist_Forest","Tropical_Grasslands",
                                                             "Savannas","Temperate_Grasslands",
                                                             "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                                             "Temperate_Mixed","Coniferous_Forests","Taiga","Tundra"))


# Plot trait distribution -------------------------------------------------

## Maybe there are some outliers, discuss this with Drew
summary(DF_trait_biome$SLA)
unique(DF_trait_biome$scrubbed_species_binomial[which(DF_trait_biome$SLA>100)])

DF_trait_biome %>% 
  group_by(Biome) %>% 
  dplyr::summarise(Median=mean(SLA,na.rm=T),
                   max=max(SLA,na.rm=T),
                   min=min(SLA,na.rm=T), sd=sd(SLA,na.rm=T)) 


cols=wes_palette("Chevalier")[c(1,4,3)]


pdf("./figs/Trait_distribution/Dist_Height.pdf",width=10)
DF_trait_biome %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Height), fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  

pdf("./figs/Trait_distribution/Dist_SLA.pdf",width=12)
DF_trait_biome %>% 
  ggplot(aes(x=SLA, y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = SLA, fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white",from=0, to=100)+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  

pdf("./figs/Trait_distribution/Dist_Seed_mass.pdf",width=12)
DF_trait_biome %>% 
  ggplot(aes(x=log(Seed_mass), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Seed_mass), fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  

pdf("./figs/Trait_distribution/Dist_Leaf_N.pdf",width=12)
DF_trait_biome %>% 
  ggplot(aes(x=Leaf_N, y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = Leaf_N, fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  



pdf("./figs/Trait_distribution/Dist_Leaf_P.pdf",width=12)
DF_trait_biome %>% 
  ggplot(aes(x=Leaf_P, y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = Leaf_P, fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  


pdf("./figs/Trait_distribution/Dist_Wood_density.pdf",width=12)
DF_trait_biome %>% 
  ggplot(aes(x=log(Wood_density), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Wood_density), fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "Woody", "No growth form info"),
                      name = "Option", guide = "legend")+
  theme_ridges(grid = FALSE)
dev.off()  



# Write data frame traits biomes ------------------------------------------
write.csv(DF_trait_biome,"./outputs/Df_traits_biomes_Sp_withRanges.csv")

