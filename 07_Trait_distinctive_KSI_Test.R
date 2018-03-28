# Libraries ---------------------------------------------------------------
library(tidyverse)

# Data --------------------------------------------------------------------
#1. Traits dataframe
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


# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass+1)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height+1)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density+1)



# KSI test ----------------------------------------------------------------
#KSI test is insensitive to data transformations
ksi_test<-function(biome,data,trait){
  
  d1<-data%>%
    filter(Biome==biome)%>%
    pull(trait)
  
  d2<-data%>%
    filter(Biome!=biome)%>%
    pull(trait)
  
  n<-c(length(d1),length(d2))
  n1n2<-sqrt(prod(n)/sum(n))
  
  fit<-suppressWarnings(ks.test(d1, d2))
  
  stat <- fit$statistic * n1n2
  
  res<-list(statistic = as.numeric(stat), p.value = fit$p.value,n = n, fit = fit)
  return(res)
}



biomes_names<-unique(Traits_Biome_Di_Ri$Biome)

Red_wides<-
Traits_Biome_Di_Ri %>% 
  filter(Ri<=0.5 & DiScale < 0.2) 

test_trait1<-lapply(biomes_names,ksi_test,Red_wides,"Leaf_P")
names(test_trait1)<-biomes_names
stat1 <- sapply(test_trait1, "[[", "statistic")
names(stat1)<-biomes_names
rev(sort(stat1))

test_trait2<-lapply(biomes_names,ksi_test,Traits_Biome_Di_Ri,"Leaf_P")
names(test_trait2)<-biomes_names
stat2 <- sapply(test_trait2, "[[", "statistic")
names(stat2)<-biomes_names
rev(sort(stat2))


