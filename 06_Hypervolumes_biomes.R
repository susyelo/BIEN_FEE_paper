# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.table("./data/base/Danilo_data/traits4susy.csv", header = TRUE)

#2. Species biomes and functional distictiveness values
sp_di_ri_biomes<-read.csv("./outputs/Biome_Di_Ri_phylo.csv")
sp_di_ri_biomes <-
  sp_di_ri_biomes %>% 
  dplyr::select(-X)

# 3. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)

# 4. 

## Merging traits and distribution dataframe ----
Traits_phylo$species<-rownames(Traits_phylo)
Biome_traits_funDi<-merge(Traits_phylo, sp_di_ri_biomes)

## Include Growth form ----
indx<-match(Biome_traits_funDi$species, Growth_form$SPECIES_STD)
Biome_traits_funDi$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]

woody<-c("Tree","Liana","Shrub","Woody epiphyte")
Biome_traits_funDi$GROWTHFORM_GEN<-ifelse(Biome_traits_funDi$GROWTHFORM_STD%in%woody,"woody","herbaceous")
Biome_traits_funDi$GROWTHFORM_GEN[which(is.na(Biome_traits_funDi$GROWTHFORM_STD))]<-NA

# Create a new levels using the grasses families

# Extract family for each species
sp_tmp<-unique(Biome_traits_funDi$species)
sp_tmp<-gsub("_"," ",sp_tmp)

family_tmp<-
  foreach(i=1:length(sp_tmp))%do% {
    print(paste("Extract",sp_tmp[i]))
    tmp<-unique(BIEN_taxonomy_species(sp_tmp[i])$scrubbed_family)
  }



# Calculate hypervolumes -------------------------------

# Scalling variables
Biome_traits_funDi$logseed_mass<-log(Biome_traits_funDi$seed_mass+1)

#Selecting and Scalling variables
scale_fun<-function(x){scale(x, center=TRUE, scale=TRUE)}

biome_name<-unique(Biome_traits_funDi$Biome)

Hyper_input_growth<-
  foreach(i=1:length(biome_name)) %do%
  {
    herbs<-Biome_traits_funDi %>% 
      filter(Biome==biome_name[i] & GROWTHFORM_GEN=="herbaceous") %>% 
      dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
      apply(2,scale_fun)
    
    woody<-Biome_traits_funDi %>% 
      filter(Biome==biome_name[i] & GROWTHFORM_GEN=="woody") %>% 
      dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
      apply(2,scale_fun)
    
    tmp<-list(woody=woody, herbs=herbs)
    
}

names(Hyper_input_growth)<-biome_name


moist_herbs<-hypervolume_box(Hyper_input$Moist_Forest$herbs, 
                       name = "Herbs")

moist_woody<-hypervolume_box(Hyper_input$Moist_Forest$woody, 
                             name = "Woody")


savanna_herbs<-hypervolume_box(Hyper_input$Savannas$herbs, 
                             name = "Herbs")

savanna_woody<-hypervolume_box(Hyper_input$Savannas$woody, 
                             name = "Woody")


plot(hypervolume_join(moist_herbs, moist_woody),num.points.max.random=6000,contour.lwd=1.5,colors=c(brewer.pal(n=2,"Set1"))
     ,show.legend=FALSE)
legend("bottomleft",legend = c("herbs","woody"),text.col=c(brewer.pal(n=3,"Set1")),bty="n",cex=1.1,text.font=2)


plot(hypervolume_join(savanna_herbs, savanna_woody),num.points.max.random=6000,contour.lwd=1.5,colors=c(brewer.pal(n=2,"Set1"))
     ,show.legend=FALSE)
legend("bottomleft",legend = c("herbs","woody"),text.col=c(brewer.pal(n=3,"Set1")),bty="n",cex=1.1,text.font=2)


savanna<-hypervolume_box(Hyper_input$Savannas, 
                       name = "Savanna")


##Divide between distinctic and redundant

FunDi_tmp<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    Biome_traits_funDi %>% 
      filter(Biome==biome_name[i]) %>% 
      mutate(Dist_redun=ifelse(FunDi>=quantile(FunDi, 0.7),"Dist","Redu"))
  }


moist_Red<-
  FunDi_tmp %>% 
  filter(Biome=="Moist_Forest" & Dist_redun=="Redu") %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  apply(2,scale_fun) %>% 
  hypervolume_box(name="Rest_moist")

moist_Dist<-
  FunDi_tmp %>% 
  filter(Biome=="Moist_Forest" & Dist_redun=="Dist") %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  apply(2,scale_fun) %>% 
  hypervolume_box(name="Dist_moist")

plot(hypervolume_join(moist_Red, moist_Dist),num.points.max.random=6000,contour.lwd=1.5,colors=c(brewer.pal(n=2,"Set1"))
     ,show.legend=FALSE)
legend("bottomleft",legend = c("Red","Dist"),text.col=c(brewer.pal(n=3,"Set1")),bty="n",cex=1.1,text.font=2)


## Using PCA
require(ggfortify)

FunDi_tmp$BiomeFu<-paste(FunDi_tmp$Biome, FunDi_tmp$Dist_redun)

FunDi_tmp %>% 
  filter(Biome=="Moist_Forest") %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  prcomp(scale=TRUE) %>% 
  autoplot(data=FunDi_tmp[which(FunDi_tmp$Biome=="Moist_Forest"),], 
           colour="Dist_redun", cex=2, frame=TRUE)
  

FunDi_tmp %>% 
  filter(Biome=="Savannas") %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  prcomp(scale=TRUE) %>% 
  autoplot(data=FunDi_tmp[which(FunDi_tmp$Biome=="Savannas"),], 
           colour="Dist_redun", cex=2, frame=TRUE)


tmp<- FunDi_tmp %>% 
  filter(Biome=="Moist_Forest" | Biome=="Dry_Forest" ) %>% 
  filter(Ri<0.5)
  
tmp %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  prcomp(scale=TRUE) %>% 
  autoplot(data=tmp, 
           colour="Biome", cex=2, frame=TRUE)



FunDi_tmp %>% 
  filter(Biome=="Taiga") %>% 
  dplyr::select(logseed_mass, wood_density, height, leaf_area) %>% 
  prcomp(scale=TRUE) %>% 
  autoplot(data=FunDi_tmp[which(FunDi_tmp$Biome=="Taiga"),], 
           colour="Dist_redun", cex=2, frame=TRUE)


## Distributions
cols=wes_palette("Chevalier")[c(1,4)]

FunDi_tmp %>% 
  ggplot(aes(x=log(height), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(height), fill = paste(Biome, Dist_redun)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Redundant", "Distinctive"),
                      name = "Functional form", guide = "legend")+
  theme_ridges(grid = FALSE)




