# Library -----------------------------------------------------------------
library(ape)
library(tidyverse)
library(picante)
library(Rphylopars)
library(ade4)
library(visdat)
library(ggtree)
library(phytools)

# Functions ---------------------------------------------------------------
source("./functions/tip_accuracy.R")

# Data --------------------------------------------------------------------
# 1. Traits data
Trait_BIEN_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv", row.names=1)
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ","_",Trait_BIEN_df$scrubbed_species_binomial)

#2. Range maps data
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv", col.names=c("Species","Y","X"))
cell_sp_biomes<-readRDS("./outputs/spPresence_biomes_all.rds")

#3. Phylogenetic data
Seed_phylo<-read.tree("./data/base/big_seed_plant_trees_v0.1/ALLMB.tre")


# Data filter and selection ----------------------------------------------------------
# Filter species that have range maps information
Trait_BIEN<-
Trait_BIEN_df %>%
  filter(scrubbed_species_binomial%in%unique(spPresence$Species)) %>%
  dplyr::select(scrubbed_species_binomial, Wood_density,Leaf_N,SLA,Seed_mass,Height,Leaf_P)

rownames(Trait_BIEN)<-Trait_BIEN$scrubbed_species_binomial

# Data exploration --------------------------------------------------------
# Number of species in the phylogeny that have some trait information
length(which(unique(Trait_BIEN$scrubbed_species_binomial)%in%Seed_phylo$tip.label))

# Trait and phylo match ----------------------------------------------------------
phylo_traits<- match.phylo.data(Seed_phylo, Trait_BIEN[,-1])
phylo_traits$data$species<-rownames(phylo_traits$data)

## Visualise the dataframe to understand the main trait gaps 
new_phylo<-ladderize(phylo_traits$phy)

indx<-match(new_phylo$tip.label,phylo_traits$data$species)
new_data<-phylo_traits$data[indx,]

tree_tmp<-ggtree(new_phylo, branch.length="none")
  
png("./supp_info/Seed_phylo_speciesWithtraits.png",width = 300)
tree_tmp
dev.off()

png("./supp_info/Missing_trait_data_phylo.png", width = 1100, height = 800)
vis_miss(new_data[,-7], sort_miss = TRUE) +
  theme(text = element_text(size=22),
        axis.text.x = element_text(size=16)) 
dev.off()

## Calculate lambda for each trait
Wood_density_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["Wood_density"]], method="lambda")
Height_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["Height"]], method="lambda")
Seed_mass_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["Seed_mass"]], method="lambda")
Leaf_N_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["Leaf_N"]], method="lambda")
Leaf_P_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["Leaf_P"]], method="lambda")
SLA_ps<-phylosig(phylo_traits$phy, phylo_traits$data[["SLA"]], method="lambda")

# Create species column and move it to the first position
phylo_traits$data<-
  phylo_traits$data %>%
  select(species, everything())


## Checking accurancy of trait imputation
# Taken from http://www.ecography.org/appendix/ecog-03480

#woody_dens<-phylo_traits$data$Wood_density
#names(woody_dens)<-phylo_traits$data$species

#accu_woody<-tip_accuracy(Tree = phylo_traits$phy, Trait = woody_dens, 
#                   method = "Rphylopars", runs = 1)



# Fill trait data using phylo info
traits_inPhylo<- phylopars(trait_data = phylo_traits$data,tree = phylo_traits$phy,
                            model = "lambda", pheno_error=FALSE,pheno_correlated = FALSE,phylo_correlated=FALSE)

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
sp_to_drop<-rownames(pca$li[which(pca$li$Axis2<(-5) & pca$li$Axis1>4),])
write.csv(sp_to_drop, "./outputs/sp_outliers.csv")

traits_completed<-
  traits_completed %>%
  filter(species%in%sp_to_drop==FALSE)
rownames(traits_completed)<-traits_completed$species

write.csv(traits_completed, "./data/processed/traits_ALLMB.csv",row.names =FALSE)



# Fill trait data using phylo info and lambda -----
traits_inPhylo2<- phylopars(trait_data = phylo_traits$data,tree = phylo_traits$phy,
                            model = "lambda", pheno_error=FALSE,pheno_correlated = FALSE,phylo_correlated=FALSE)

traits_completed<-as.data.frame(traits_inPhylo2$anc_recon[1:length(phylo_traits$phy$tip.label),])

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
sp_to_drop<-rownames(pca$li[which(pca$li$Axis2<(-5) & pca$li$Axis1>4),])
write.csv(sp_to_drop, "./outputs/sp_outliers.csv")

traits_completed<-
  traits_completed %>%
  filter(species%in%sp_to_drop==FALSE)
rownames(traits_completed)<-traits_completed$species

write.csv(traits_completed, "./data/processed/traits_ALLMB_lambda.csv",row.names =FALSE)
