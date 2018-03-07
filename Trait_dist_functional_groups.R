# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)


# data --------------------------------------------------------------------
DF_trait_biome<-read.csv("./outputs/Df_traits_biomes_Sp_withRanges.csv")


DF_trait_biome_sub<-
  DF_trait_biome%>%
  select(scrubbed_species_binomial,SLA,Seed_mass,Height,Leaf_P,Leaf_N,Wood_density) %>% 
  distinct

rownames(DF_trait_biome_sub)<-DF_trait_biome_sub$scrubbed_species_binomial

# Compute distance matrix of trait between each pair of species 
Dist_matrix<-compute_dist_matrix(DF_trait_biome_sub[,-1],metric="euclidean",center = TRUE,
                                 scale = TRUE) ## This can take a while


## Check this file for further analysis
# /Users/echeverrialondono1/Git_Repos/BIEN_project/Functional_Div/4__Dominants_Subordinates_Dist.R
