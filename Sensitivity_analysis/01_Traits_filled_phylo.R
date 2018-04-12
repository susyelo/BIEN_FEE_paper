# Library -----------------------------------------------------------------
library(ape)
library(tidyverse)
library(picante)
library(Rphylopars)
library(ade4)
library(visdat)
library(raster)
library(letsR)


# Functions ---------------------------------------------------------------
source("./functions/tip_accuracy.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")

# Data --------------------------------------------------------------------
# 1. Traits data
Trait_BIEN_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv", row.names=1)
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ","_",Trait_BIEN_df$scrubbed_species_binomial)

#2. Range maps data
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv", col.names=c("Species","Y","X"))

#3. Phylogenetic data
Seed_phylo<-read.tree("./data/base/big_seed_plant_trees_v0.1/ALLMB.tre")


# 4. Richness raster
r<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

# Data manipulation ----------------------------------------------------------
Trait_BIEN_df<-
Trait_BIEN_df %>% 
  dplyr::select(scrubbed_species_binomial,Height,Seed_mass,Wood_density,SLA,Leaf_N,Leaf_P)

rownames(Trait_BIEN_df)<-Trait_BIEN_df$scrubbed_species_binomial

# Trait and phylo match ----------------------------------------------------------
phylo_traits<- match.phylo.data(Seed_phylo, Trait_BIEN_df[,-1])
phylo_traits$data$species<-rownames(phylo_traits$data)

# Create species column and move it to the first position
phylo_traits$data<-
  phylo_traits$data %>%
  select(species, everything())

# Fill trait data by default -----------------------------------------------------------
traits_inPhylo<- phylopars(trait_data = phylo_traits$data,tree = phylo_traits$phy,
                     pheno_error = F, phylo_correlated = F, pheno_correlated = F)

traits_completed<-as.data.frame(traits_inPhylo$anc_recon[1:length(phylo_traits$phy$tip.label),])

# Drop outliers 
pca <- dudi.pca(traits_completed,
                scannf = F, nf = 5)

plot(pca$li[,1:2])
pca$li[which(pca$li$Axis1<(-5)),]

# Dropping Cocos_nucifera
traits_completed$species<-rownames(traits_completed)

traits_completed<-
  traits_completed %>%
  filter(species!="Cocos_nucifera") %>%
  droplevels()


# Testing again for outliers
rownames(traits_completed)<-traits_completed$species
pca <- dudi.pca(traits_completed[,-7],
                scannf = F, nf = 5)

plot(pca$li[,1:2])

## Extreme axis values based mainly on extreme seed mass values
sp_to_drop<-rownames(pca$li[which(pca$li$Axis1<(-4) & pca$li$Axis2>10),])

traits_completed<-
  traits_completed %>%
  filter(species%in%sp_to_drop==FALSE)

# Maps with default input -------------------------------------------------
## Species with range maps
rownames(traits_completed)<-traits_completed$species

traits_completed_map<-
  traits_completed %>%
  filter(species%in%unique(spPresence$Species))


## Richness map
spRichness = splistToRichness(spPresence,traits_completed_map$species)
spRichness[spRichness==0]<-NA

Trait_richness_plot<-spplot(spRichness, 
                            main=paste("Species with traits",
                                       length(traits_completed_map$species),"sp"))

## presence/absence matrix
spMatrix<-splistToMatrix(spPresence,traits_completed_map$species)

# letsR presenceAbsence objects
pb_spp<-list(Presence_and_Absence_Matrix = spMatrix, 
             Richness_Raster = spRichness, Species_name = colnames(spMatrix))
class(pb_spp) <- "PresenceAbsence"


# Height map --------------------------------------------------------------
o<-match(pb_spp$Species_name,traits_completed_map$species)
height<-traits_completed_map$Height[o]

height_ras<-lets.maplizer.mod(pb_spp,height,pb_spp$Species_name,func = mean, ras=TRUE)
height_plot<-spplot(log10(height_ras$Raster))


# Maps with lambda method input -------------------------------------------

