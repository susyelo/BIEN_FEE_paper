# Library -----------------------------------------------------------------
library(ape)
library(tidyverse)
library(picante)
library(Rphylopars)
library(ade4)
library(visdat)
library(ggtree)
library(phytools)
library(RColorBrewer)
library(foreach)

# Functions ---------------------------------------------------------------
source("./functions/tip_accuracy.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")

# Data --------------------------------------------------------------------
# 1. Traits data
Trait_BIEN_df<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv", row.names=1)
Trait_BIEN_df$scrubbed_species_binomial<-gsub(" ","_",Trait_BIEN_df$scrubbed_species_binomial)

#2. Range maps data
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv",header = FALSE, col.names=c("Species","Y","X"))
cell_sp_biomes<-readRDS("./outputs/spPresence_biomes_all.rds")

#3. Phylogenetic data
Seed_phylo<-read.tree("./data/base/big_seed_plant_trees_v0.1/ALLMB.tre")

#4. Total_richness raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

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


# Geographic distribution of traits ---------------------------------------
## Richness map and among biomes
p<-cell_sp_biomes %>% 
  group_by(cells) %>% 
  summarise(Richness=n_distinct(Species), biomes=unique(biomes)) %>% 
  ggplot(aes(x=reorder(biomes, Richness, FUN=mean), y=Richness)) +
  geom_boxplot() 

p +  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Biomes") + ylab("Richness")

## Presence/absence matrix of species
## Total richness
spRichness = splistToRichness(spPresence,unique(cell_sp_biomes$Species))
spRichness[spRichness==0]<-NA

Richness_plot<-spplot(spRichness)

## Richness of species with at least one trait
spRichness_traits = splistToRichness(spPresence,unique(Trait_BIEN$scrubbed_species_binomial))
spRichness_traits[spRichness_traits==0]<-NA

Richness_traits_plot<-spplot(spRichness_traits)

## proportion of species with traits
prop_richnnes<-spRichness_traits/spRichness

my.palette=c('#ffffcc','#c2e699','#78c679','#31a354','#006837')
spl_tmp<- list("sp.lines", as(biome_shp, "SpatialLines"), col=alpha("dimgrey",0.4))


png("./supp_info/Prop_species_traits.png")
spplot(prop_richnnes,col.regions = my.palette, at = seq(0, 0.5, 0.1), 
       par.settings = list(fontsize = list(text = 40)),sp.layout=spl_tmp)
dev.off()


## proportion per traits

trait_names<-names(Trait_BIEN)[-1]

prop_Sp_trait_rasters<-foreach(i=1:length(trait_names))%do%{
  
  # 1. Extract list of species with available trait info
  sp_trait<-
    Trait_BIEN %>% 
    dplyr::select(scrubbed_species_binomial,mytrait=trait_names[i]) %>% 
    filter(!is.na(mytrait)) %>% 
    pull(scrubbed_species_binomial)
  
  # 2. Estimate richness per grid cell of species with trait info
  spRichness_trait = splistToRichness(spPresence,sp_trait)
  spRichness_trait[spRichness_trait==0]<-NA
  
  # 3. Raster with the proportion of species with sampled trait
  prop_richnnes_trait<-spRichness_trait/spRichness
  
  prop_richnnes_trait
}

names(prop_Sp_trait_rasters)<-trait_names


p1<-spplot(prop_Sp_trait_rasters$Leaf_P,col.regions = my.palette,  at = seq(0, 0.5, 0.1),
       par.settings = list(fontsize = list(text = 35)),sp.layout=spl_tmp,colorkey=FALSE)

p2<-spplot(prop_Sp_trait_rasters$Leaf_N,col.regions = my.palette, at = seq(0, 0.5, 0.1),
       par.settings = list(fontsize = list(text = 35)),sp.layout=spl_tmp,colorkey=FALSE)

p3<-spplot(prop_Sp_trait_rasters$Wood_density,col.regions = my.palette, at = seq(0, 0.5, 0.1),
       par.settings = list(fontsize = list(text = 40)),sp.layout=spl_tmp,colorkey=FALSE)

p4<-spplot(prop_Sp_trait_rasters$SLA,col.regions = my.palette,at = seq(0, 0.5, 0.1),
       par.settings = list(fontsize = list(text = 35)),sp.layout=spl_tmp,colorkey=FALSE)

p5<-spplot(prop_Sp_trait_rasters$Height,col.regions = my.palette, cuts=5, at = seq(0, 0.5, 0.1),
       par.settings = list(fontsize = list(text = 35)),sp.layout=spl_tmp,colorkey=FALSE)


p6 <- c(p1,p2,p3,p4,p5, layout=c(6,1))

png("./supp_info/all_prop_traits.png", width = 1000, height = 400)
p6
dev.off()

pdf("./supp_info/all_prop_traits.pdf", width = 12, height = 5)
p6
dev.off()

pdf("./supp_info/Prop_species_traits_Seed_mass.pdf")
spplot(prop_Sp_trait_rasters$Seed_mass,col.regions = my.palette, at = seq(0, 0.5, 0.1),
           par.settings = list(fontsize = list(text = 35)),sp.layout=spl_tmp)
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
