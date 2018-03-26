# libraries ---------------------------------------------------------------
library(tidyverse)
library(BIEN)

# Functions ---------------------------------------------------------------
source("./functions/Plot_Half_Pies.R")

# data --------------------------------------------------------------------
# 1. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

#2. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)

# 3. Trait values for species
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

# 4. Merging Distinctiveness dataframes and traits dataframe
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)


# Include growth form -----------------------------------------------------
indx<-match(Traits_Biome_Di_Ri$species,Growth_form$SPECIES_STD)
Traits_Biome_Di_Ri$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]

# Include general forms
woody<-c("Tree","Liana","Shrub","Woody epiphyte")
Traits_Biome_Di_Ri$GROWTHFORM_GEN<-ifelse(Traits_Biome_Di_Ri$GROWTHFORM_STD%in%woody,"woody","herbaceous")
Traits_Biome_Di_Ri$GROWTHFORM_GEN[which(is.na(Traits_Biome_Di_Ri$GROWTHFORM_STD))]<-NA

# Include grasses as a Growth form


# Extract dominant growth forms per biome ---------------------------------
biome_name<-unique(Traits_Biome_Di_Ri$Biome)

Total<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Traits_Biome_Di_Ri %>% 
      filter(Biome==biome_name[i]) %>% 
      group_by(GROWTHFORM_STD) %>% 
      dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
      mutate(Dist="Total_prop",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,1)) %>% 
      filter(prop > 1)
  }

#Using the ScaleUi values produced the same results as the ScaleDi
Redundant_widespread<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Traits_Biome_Di_Ri %>% 
      filter(Biome==biome_name[i]) %>% 
      filter(Ri<=0.5 & DiScale < 0.2) %>% 
      group_by(GROWTHFORM_STD) %>% 
      dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
      mutate(Dist="Redun_wides",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,1)) %>% 
      filter(prop > 2.1)
  }

# The number of species per biome in this category is not significant, therefore I decided not included in the analysis

Distinctive_widespread<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Traits_Biome_Di_Ri %>% 
      filter(Biome==biome_name[i]) %>% 
      filter(Ri<=0.5 & DiScale > 0.1) %>% 
      group_by(GROWTHFORM_STD) %>% 
      dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
      mutate(Dist="Dist_wides",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,1)) %>% 
      filter(prop > 2.1)
  }


dir.create("./figs/Growth_forms")

png("./figs/Growth_forms/Total_proportion.png",width = 500)

ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = Total,
                    stat="identity")+
  coord_flip() +
  ggtitle("Total_proportion")

dev.off()


png("./figs/Growth_forms/Redundant_widespread.png",width = 500)

ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = Redundant_widespread,
                    stat="identity")+
  coord_flip() +
  ggtitle("Redundant widespread")

dev.off()


png("./figs/Growth_forms/Distinctive_widespread.png",width = 500)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = Distinctive_widespread,
                    stat="identity")+
  coord_flip() +
  ggtitle("Distinctive widespread")

dev.off()


## Semicircles pies

tmp<-
Total %>%
  filter(Biome=="Moist_Forest")
cols=wes_palette("Darjeeling",n_distinct(tmp$GROWTHFORM_STD))

library(wesanderson)  
moist<-parlDiag(tmp$GROWTHFORM_STD, tmp$N_sp, cols =cols, repr = "absolute")
moist2<-parlDiag(tmp$GROWTHFORM_STD, tmp$N_sp, cols =cols, repr = "absolute")


library(cowplot)
plot_grid(moist, moist2, labels=c("A", "B"), ncol = 2, nrow = 1)

# Trait distribution: Redundant and widespread species --------------------

## Total distribution
Traits_Biome_Di_Ri %>%
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges2(aes(x = log(Height), fill = Biome),calc_ecdf = TRUE,
                       rel_min_height = 0.01,scale = 1)+theme_ridges()

Traits_Biome_Di_Ri %>% 
  filter(Ri<=0.5 & DiScale > 0.1) %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges2(aes(x = log(Height), fill = Biome),calc_ecdf = TRUE,
                       rel_min_height = 0.01,scale = 1)+theme_ridges()+ 
  theme_minimal(base_size = 14) + theme(axis.text.y = element_text(vjust = 0)) 




Traits_Biome_Di_Ri %>% 
  filter(Ri<=0.5 & DiScale > 0.1) %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Height), fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "No information", "Woody"),
                      name = "Growth form", guide = "legend")+
  theme_ridges(grid = FALSE)
