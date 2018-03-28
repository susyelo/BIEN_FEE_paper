# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)
library(factoextra)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

# 3. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)

Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome, 
                                 levels=rev(c("Moist_Forest","Savannas","Tropical_Grasslands",
                                              "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                              "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                              "Taiga","Tundra")))

# Calculate hypervolumes -------------------------------

# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass+1)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height+1)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density+1)

#Selecting and Scalling variables
Traits_Biome_Di_Ri<-
Traits_Biome_Di_Ri %>% 
  group_by(Biome) %>% 
  mutate(Scaled_logWood_density=scale(logWoodDensity),
         Scaled_Leaf_N=scale(Leaf_N),
         Scaled_SLA=scale(SLA),
         Scaled_Leaf_P=scale(Leaf_P),
         Scaled_logSeed_mass=scale(logseed_mass),
         Scaled_logHeight=scale(logHeight))
  

moist<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Biome=="Moist_Forest") %>% 
  dplyr::filter(Ri<=0.5 & DiScale < 0.2) %>% 
  select(contains("Scaled"))

moist_hb<-hypervolume_box(moist[,-1],name = "Moist_Forest")

Dry<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Biome=="Dry_Forest") %>% 
  dplyr::filter(Ri<=0.5 & DiScale < 0.2) %>% 
  select(contains("Scaled"))

Dry_hb<-hypervolume_box(Dry[,-1],name = "Dry_Forest")


plot(hypervolume_join(moist_hb, Dry_hb),contour.lwd=1.5,colors=c(brewer.pal(n=3,"Set1"))
     ,show.legend=TRUE)


sim1<-
hypervolume_overlap_statistics(
  
  hypervolume_set(
    moist_hb,Dry_hb,
    check.memory = FALSE
  )
)

## With the total values
moist_tot<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Biome=="Moist_Forest") %>%
  select(contains("Scaled"))

moist_hb_tot<-hypervolume_box(moist_tot[,-1],name = "Moist_Forest")

Dry_tot<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Biome=="Dry_Forest") %>% 
  select(contains("Scaled"))

Dry_hb_tot<-hypervolume_box(Dry_tot[,-1],name = "Dry_Forest")


plot(hypervolume_join(moist_hb_tot, Dry_hb_tot),contour.lwd=1.5,colors=c(brewer.pal(n=3,"Set1"))
     ,show.legend=TRUE)
z_nolog <- hypervolume(data_twoboxes,bandwidth=estimate_bandwidth(data_twoboxes),quantile=0.5,reps=1000)

sim2<-
  hypervolume_overlap_statistics(
    hypervolume_set(
      moist_hb_tot,Dry_hb_tot,
      check.memory = FALSE
    )
  )

## Test pca

traits_tot<-
  Traits_Biome_Di_Ri %>% 
  select(Leaf_N,Leaf_P,SLA,logseed_mass,logHeight,logWoodDensity)

trait_pca <- prcomp(traits_tot[,-1],scale. = TRUE)
summary(trait_pca)

fviz_pca_var(trait_pca)




dat_sub<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Ri<=0.5 & DiScale < 0.2)

traits_sub<-
  dat_sub %>% 
  select(Leaf_N,Leaf_P,SLA,logseed_mass,logHeight,logWoodDensity)

trait_pca_sub <- prcomp(traits_sub[,-1],scale. = TRUE)
summary(trait_pca_sub)

fviz_pca_var(trait_pca_sub)

fviz_pca_var(trait_pca_sub, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=10) +
  theme_minimal()

fviz_screeplot(trait_pca_sub, addlabels = TRUE, ylim = c(0, 50))

iris.pca <- PCA(traits_sub[,-1], graph = FALSE)
fviz_pca_ind(iris.pca, label="none", habillage=dat_sub$Biome,
             addEllipses=TRUE, ellipse.level=0.95)
