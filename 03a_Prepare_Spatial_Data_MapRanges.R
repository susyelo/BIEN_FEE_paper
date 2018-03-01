# Libraries ---------------------------------------------------------------
library(letsR)
library(raster)
library(gridExtra)
library(tidyverse)


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


# Total richness ----------------------------------------------------------
### Total richness for species with ranges maps in BIEN 2.0
plot(r)
r[r==0]<-NA
Total_richness_plot<-spplot(r, main=paste("Total richness:",n_distinct(d$Species),"sp"))


# Richness for species with trait values ----------------------------------

# process trait dataframe
trait_df$scrubbed_species_binomial<-gsub(" ","_",trait_df$scrubbed_species_binomial)

trait_df_wide<-
  trait_df  %>% 
  spread(trait_name,trait_value) %>% 
  filter(scrubbed_species_binomial%in%d$Species)


# matrix and richness raster
TraitSpecies = unique(trait_df_wide$scrubbed_species_binomial)

## It include NAs for the species that do not have range maps available
spMatrix = splistToMatrix(d,TraitSpecies)

# Richness using only the species with trait values -----------------------
spRichness = splistToRichness(d,TraitSpecies)
spRichness[spRichness==0]<-NA

Trait_richness_plot<-spplot(spRichness, main=paste("Species with traits",length(TraitSpecies),"sp"))


## Richness for woody vs herbaceous species
## Woody species
woody_sp<-trait_df_wide %>% 
  filter(GROWTHFORM_GEN=="woody") %>% 
  dplyr::select(scrubbed_species_binomial) %>% 
  unique

Woody_Richness = splistToRichness(d,woody_sp$scrubbed_species_binomial)
Woody_Richness[Woody_Richness==0]<-NA
Woody_plot<-spplot(Woody_Richness, main=paste("Woody",n_distinct(woody_sp$scrubbed_species_binomial),"sp"))


## Herbaceous species
herbaceous_sp<-trait_df_wide %>% 
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


