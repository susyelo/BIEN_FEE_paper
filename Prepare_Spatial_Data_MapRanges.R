# Libraries ---------------------------------------------------------------
library(letsR)
library(raster)
library(gridExtra)
library(tidyverse)
library(foreach)
library(maptools)


# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")


# Data --------------------------------------------------------------------
# Trait data
trait_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv")

# Range maps
d<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv")
names(d) = c("Species","Y","X")

# Richness raster
r<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

# Biomes polygons
load("./data/base/Danilo_data/biomes_shp4susy.RData")


# Total richness ----------------------------------------------------------
### Total richness for species with ranges maps in BIEN 2.0
plot(r)
r[r==0]<-NA
Total_richness_plot<-spplot(r, main=paste("Total richness:",n_distinct(d$Species),"sp"))


# Richness for species with trait values ----------------------------------

# process trait dataframe
trait_df$scrubbed_species_binomial<-gsub(" ","_",trait_df$scrubbed_species_binomial)

trait_df_range<-
  trait_df  %>% 
  filter(scrubbed_species_binomial%in%d$Species)

# matrix and richness raster
TraitSpecies = unique(trait_df_range$scrubbed_species_binomial)

## It include NAs for the species that do not have range maps available
spMatrix = splistToMatrix(d,TraitSpecies)

# Richness using only the species with trait values -----------------------
spRichness = splistToRichness(d,TraitSpecies)
spRichness[spRichness==0]<-NA

Trait_richness_plot<-spplot(spRichness, main=paste("Species with traits",length(TraitSpecies),"sp"))


## Richness for woody vs herbaceous species
## Woody species
woody_sp<-trait_df_range %>% 
  filter(GROWTHFORM_GEN=="woody") %>% 
  dplyr::select(scrubbed_species_binomial) %>% 
  unique

Woody_Richness = splistToRichness(d,woody_sp$scrubbed_species_binomial)
Woody_Richness[Woody_Richness==0]<-NA
Woody_plot<-spplot(Woody_Richness, main=paste("Woody",n_distinct(woody_sp$scrubbed_species_binomial),"sp"))


## Herbaceous species
herbaceous_sp<-trait_df_range %>% 
  filter(GROWTHFORM_GEN=="herbaceous") %>% 
  dplyr::select(scrubbed_species_binomial) %>% 
  unique

herbaceus_Richness = splistToRichness(d,herbaceous_sp$scrubbed_species_binomial)
herbaceus_Richness[herbaceus_Richness==0]<-NA
herbaceus_plot<-spplot(herbaceus_Richness, main=paste("Herbaceus",n_distinct(herbaceous_sp$scrubbed_species_binomial),"sp"))


pdf("./figs/Rangemaps_richness.pdf")
grid.arrange(Total_richness_plot,Trait_richness_plot,
             Woody_plot,herbaceus_plot,
             ncol=2,
             nrow=2)
dev.off()


# Writing data ------------------------------------------------------------
writeRaster(spRichness,"./data/processed/TraitRichness_raster.tif")
save(spMatrix, file="./data/processed/PresAbs_matrix_TraitRange.RData")


## Species richness per biome

biome_poly<-list(coniferous,dry,grasslands,
                 moist,montane,prairies,savanna,taiga,
                 temperate.mixed,tropical.mixed,tundra)

names(biome_poly)<-c("coniferous","dry","grasslands",
                     "moist","montane","prairies","savanna","taiga",
                     "temperate.mixed","tropical.mixed","tundra")

biome_richness<-foreach(i=1:length(biome_poly))  %do% {
  
  unlist(raster::extract(r,biome_poly[[i]]))

}

names(biome_richness)<-names(biome_poly)


## transform data
require(reshape2)
biome_richness_df <- melt(biome_richness)
colnames(biome_richness_df)<-c("N_Species","Biome")
biome_richness_df<-na.omit(biome_richness_df)

## Plotting Violin plots for each biome
p <- ggplot(biome_richness_df, 
            aes(x=reorder(Biome, N_Species, FUN=mean), y=N_Species)) +
  geom_violin(trim=TRUE)

p + stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Biomes") + ylab("Richness")


### Species richness with traits per biome ----

Biome_richness_trait<-foreach(i=1:length(biome_poly))  %do% {
  
  unlist(raster::extract(spRichness,biome_poly[[i]]))
  
}

names(Biome_richness_trait)<-names(biome_poly)


## transform data
require(reshape2)
biome_richness_trait_df <- melt(Biome_richness_trait)
colnames(biome_richness_trait_df)<-c("N_Species","Biome")
biome_richness_trait_df<-na.omit(biome_richness_trait_df)

## Plotting Violin plots for each biome
p <- ggplot(biome_richness_trait_df, 
            aes(x=reorder(Biome, N_Species, FUN=mean), y=N_Species)) +
  geom_violin(trim=TRUE)

p + stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Biomes") + ylab("Richness")

