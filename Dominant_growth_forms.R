# libraries ---------------------------------------------------------------
library(tidyverse)

# data --------------------------------------------------------------------
# 1. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

#2. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)


# Extract dominant growth forms per biome ---------------------------------



dist<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Biome_Di_Ri %>% 
      filter(Biome==biome_name[i]) %>% 
      filter(FunDi>=quantile(FunDi, 0.7)) %>% 
      group_by(GROWTHFORM_STD) %>% 
      dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
      mutate(Dist="Dist",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,2))
  }