# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)
library(factoextra)
library(dendextend)
library(wesanderson)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/Biomes_hypervolumes_fun.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

#3. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)


# 4. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)


## Rename and order biomes
Traits_Biome_Di_Ri$Biome<-recode(Traits_Biome_Di_Ri$Biome,Moist_Forest="Moist",
                                     Savannas="Savannas",
                                     Tropical_Grasslands="Trop_Grass",
                                     Dry_Forest="Dry",
                                     Xeric_Woodlands="Xeric",
                                     Mediterranean_Woodlands="Mediterranean",
                                     Temperate_Grasslands="Temp_Grass",
                                     Temperate_Mixed="Temp_Mixed",
                                     Coniferous_Forests="Coniferous",
                                     Taiga="Taiga",
                                     Tundra="Tundra")

Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome, 
                                 levels=c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra"))

# Calculate hypervolumes -------------------------------

# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density)
Traits_Biome_Di_Ri$sqrtSLA<-sqrt(Traits_Biome_Di_Ri$SLA)

#Selecting and Scalling variables
Traits_Biome_Di_Ri<-
Traits_Biome_Di_Ri %>% 
  group_by(Biome) %>% 
  mutate(Scaled_logWood_density=scale(logWoodDensity),
         Scaled_Leaf_N=scale(Leaf_N),
         Scaled_SLA=scale(sqrtSLA),
         Scaled_Leaf_P=scale(Leaf_P),
         Scaled_logSeed_mass=scale(logseed_mass),
         Scaled_logHeight=scale(logHeight))
  

# Redundant and widespread species hypervolumes ---------------------------
## Ordering traits
biome_names<-levels(Traits_Biome_Di_Ri$Biome)

## Hypervolumes for all species
Total_hypervol<-
  Traits_Biome_Di_Ri %>% 
  dplyr::select(Biome,contains("Scaled")) %>% 
  Biomes_hypervolume(., biome_names)

saveRDS(Total_hypervol, "./outputs/Total_hypervolumes.rds")

## Hypervolumes for widespread and redundant species 
Redun_Wides_hypervol<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Ri<=0.5 & DiScale < 0.2) %>% 
  dplyr::select(Biome,contains("Scaled")) %>% 
  Biomes_hypervolume(biome_names)

saveRDS(Redun_Wides_hypervol, "./outputs/ReduntWides_hypervolumes.rds")


## Plot hypervolumes per category
png("./figs/hypervolumes_clusters/Total_Moist_Dry_Savanna.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Total_hypervol$Moist_Forest, 
    Total_hypervol$Dry_Forest,
    Total_hypervol$Savannas
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/Total_Moist_Temperated_Mixed_Conifers.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Total_hypervol$Moist_Forest, 
    Total_hypervol$Temperate_Mixed,
    Total_hypervol$Coniferous_Forests
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/RedWides_Moist_Dry_Savanna.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist_Forest, 
    Redun_Wides_hypervol$Dry_Forest,
    Redun_Wides_hypervol$Savannas
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/RedWides_Moist_Temperated_Mixed_Conifers.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist_Forest, 
    Redun_Wides_hypervol$Temperate_Mixed,
    Redun_Wides_hypervol$Coniferous_Forests
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()


## Calculate Hypervolume similarity using Sorense's index -----
## With Redundant species
redun_Sim<-similarity_hypervol(Redun_Wides_hypervol)
fit_red <-hclust(as.dist(1-redun_Sim))


dend_red<-
  fit_red %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti1")[3]) %>% 
  set("branches_lwd", 4) %>% 
  set("labels_cex", 1.5)


## With the total species
Total_Sim<-similarity_hypervol(Total_hypervol)
fit_total <-hclust(as.dist(1-Total_Sim))

labels(fit_total)<-c("Trop_Grass","Moist","Savannas","Dry","Taiga","Tundra","Xeric",
                     "Mediterranean","Temp_Grass","Temp Mixed","Coniferous")

dend_total<-
  fit_total %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti1")[1]) %>% 
  set("branches_lwd", 4) %>% 
  set("labels_cex", 1.5)


#dir.create("./figs/hypervolumes_clusters")
pdf("./figs/hypervolumes_clusters/Redundant_Sorensen.pdf", height = 9.4, width = 9.1)
circlize_dendrogram(dend_red,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

pdf("./figs/hypervolumes_clusters/Total_Sorensen.pdf",height = 9.4, width = 9.1)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

# Hypervolumes for climatic categories ------------------------------------
tropical<-c("Moist_Forest","Dry_Forest","Tropical_Grasslands","Savannas")
temperate<-c("Temperate_Grasslands","Coniferous_Forests","Temperate_Mixed","Mediterranean_Woodlands")
cold<-c("Taiga","Tundra")

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%tropical,"Tropical",NA)
Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%temperate,"Temperate",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%cold,"Cold",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes[which(Traits_Biome_Di_Ri$Biome=="Xeric_Woodlands")]<-"Xeric"


## The function needs a variable called "Biome" 
## so I created a temporal dataframe with the Biome variables using the climate classification
Traits_Biome_Di_Ri_tmp<-Traits_Biome_Di_Ri
Traits_Biome_Di_Ri_tmp$Biome<-Traits_Biome_Di_Ri_tmp$ClimBiomes

biome_names<-unique(Traits_Biome_Di_Ri_tmp$Biome)


Climate_Biome_Di_Ri<-
  Traits_Biome_Di_Ri_tmp %>% 
  dplyr::select(contains("Scaled"))

Climatic_hypervol<-Biomes_hypervolume(Climate_Biome_Di_Ri,biome_names)
saveRDS(Climatic_hypervol,"./outputs/Climatic_hypervolumes.rds")


Climatic_Sim<-similarity_hypervol(Climatic_hypervol)
fit_Climatic <-hclust(as.dist(1-Climatic_Sim))

dend_climatic<-
  fit_Climatic %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti1")[1]) %>% 
  set("branches_lwd", 4)



png("./figs/hypervolumes_clusters/Climatic_hypervolumes_total.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Climatic_hypervol$Tropical, 
    Climatic_hypervol$Xeric,
    Climatic_hypervol$Temperate,  
    Climatic_hypervol$Cold
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  cex.data=1,cex.axis=1,cex.names=1,cex.legend=2,
  show.random=FALSE
)
dev.off()

### Hypervolumes for Dominant growth forms




