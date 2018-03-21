# Library -----------------------------------------------------------------
library(ape)
library(tidyverse)
library(picante)
library(Rphylopars)
library(ade4)

# Data --------------------------------------------------------------------
# 1. Traits data
Trait_BIEN_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv", row.names=1)
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ","_",Trait_BIEN_df$scrubbed_species_binomial)

#2. Phylogenetic data
Seed_phylo<-read.tree("./data/base/big_seed_plant_trees_v0.1/ALLMB.tre")

# Number of species in the phylogeny that have some trait information
length(which(Seed_phylo$tip.label%in%Trait_BIEN_df$scrubbed_species_binomial))


# Data selection ----------------------------------------------------------
Trait_BIEN<-
  Trait_BIEN_df %>% 
  dplyr::select(Wood_density,Leaf_N,SLA,Seed_mass,Height,Leaf_P)

rownames(Trait_BIEN)<-Trait_BIEN_df$scrubbed_species_binomial

# Trait and phylo match ----------------------------------------------------------
phylo_traits<- match.phylo.data(Seed_phylo, Trait_BIEN)

# Create species column and move it to the first position
phylo_traits$data$species<-rownames(phylo_traits$data)

phylo_traits$data<-
  phylo_traits$data %>% 
  select(species, everything())

traits_inPhylo<- phylopars(trait_data = phylo_traits$data,tree = phylo_traits$phy,
                     pheno_error = F, phylo_correlated = F, pheno_correlated = F)


traits_completed<-as.data.frame(traits_inPhylo$anc_recon[1:length(phylo_traits$phy$tip.label),])

# Drop outliers -----------------------------------------------------------
pca <- dudi.pca(traits_completed,
                scannf = F, nf = 5)

plot(pca$li[,1:2])
pca$li[which(pca$li$Axis1>5),]

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
sp_to_drop<-rownames(pca$li[which(pca$li$Axis2>10 & pca$li$Axis1>5),])
write.csv(sp_to_drop, "./outputs/sp_outliers.csv")

traits_completed<-
  traits_completed %>% 
  filter(species%in%sp_to_drop==FALSE)
rownames(traits_completed)<-traits_completed$species

write.csv(traits_completed, "./data/processed/traits_ALLMB.csv",row.names =FALSE)
