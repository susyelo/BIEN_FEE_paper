# libraries ---------------------------------------------------------------
library(tidyverse)
library(BIEN)
library(foreach)
library(wesanderson)
library(ggridges)

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

# 5. Grasses_information
family_info<-read.csv("./data/base/family_completed.csv", row.names = 1)
family_info$species<-gsub(" ","_",family_info$species)
grass_families<-read.csv("./data/base/Grass_and_GrassLike_families.csv")

# Include growth form -----------------------------------------------------
indx<-match(Traits_Biome_Di_Ri$species,Growth_form$SPECIES_STD)
Traits_Biome_Di_Ri$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]

# Include general forms
woody<-c("Tree","Liana","Shrub","Woody epiphyte")
Traits_Biome_Di_Ri$GROWTHFORM_GEN<-ifelse(Traits_Biome_Di_Ri$GROWTHFORM_STD%in%woody,"woody","herbaceous")
Traits_Biome_Di_Ri$GROWTHFORM_GEN[which(is.na(Traits_Biome_Di_Ri$GROWTHFORM_STD))]<-NA

# Include grasses as a Growth form
indx<-match(Traits_Biome_Di_Ri$species,family_info$species)
Traits_Biome_Di_Ri$Family<-family_info$family[indx]

grasses<-unique(grass_families$Family)
Traits_Biome_Di_Ri$GROWTHFORM_STD<-as.character(Traits_Biome_Di_Ri$GROWTHFORM_STD)
Traits_Biome_Di_Ri$GROWTHFORM_STD[which(Traits_Biome_Di_Ri$Family%in%grasses)]<-"Grass"

# Include the following into the "hebs" category
joinGF<-c("Vine","Liana","Non-woody epiphyte","Aquatic","Parasite")
Traits_Biome_Di_Ri$GROWTHFORM_STD[which(Traits_Biome_Di_Ri$GROWTHFORM_STD%in%joinGF)]<-"Herb"
Traits_Biome_Di_Ri$GROWTHFORM_STD<-as.factor(Traits_Biome_Di_Ri$GROWTHFORM_STD)


# Change biome levels orders 
Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome, 
                                 levels=rev(c("Moist_Forest","Savannas","Tropical_Grasslands",
                                          "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                          "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                          "Taiga","Tundra")))

# Change Growth forms levels orders
Traits_Biome_Di_Ri$GROWTHFORM_GEN[is.na(Traits_Biome_Di_Ri$GROWTHFORM_GEN)]<-"Unknown"

Traits_Biome_Di_Ri$GROWTHFORM_STD<-as.character(Traits_Biome_Di_Ri$GROWTHFORM_STD)
Traits_Biome_Di_Ri$GROWTHFORM_STD[is.na(Traits_Biome_Di_Ri$GROWTHFORM_STD)]<-"Unknown"

Traits_Biome_Di_Ri$GROWTHFORM_GEN<-factor(Traits_Biome_Di_Ri$GROWTHFORM_GEN, 
                                          levels=c("woody","herbaceous","Unknown"))


Traits_Biome_Di_Ri$GROWTHFORM_STD<-factor(Traits_Biome_Di_Ri$GROWTHFORM_STD,
                                          levels=c("Tree","Shrub","Herb",
                                                   "Grass","Unknown"))


# Correcting growth form in Tundra ----------------------------------------
sp_tundra_nonWoody<-
Traits_Biome_Di_Ri %>% 
  filter(Biome=="Tundra"&log(Height)<(-2)&GROWTHFORM_GEN=="woody") %>% 
  select(species)

Traits_Biome_Di_Ri$GROWTHFORM_GEN[which(Traits_Biome_Di_Ri$species%in%sp_tundra_nonWoody$species)]<-"herbaceous"
Traits_Biome_Di_Ri$GROWTHFORM_STD[which(Traits_Biome_Di_Ri$species%in%sp_tundra_nonWoody$species)]<-"Herb"


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


col_GF<-c(wes_palette("Cavalcanti")[c(2:4,1)],"grey")
 
dir.create("./figs/Growth_forms")

### Both graphics into one
Total$Tmnt<-"Total"
Redundant_widespread$Tmnt<-"RedWid"

new_df<-rbind(Total,Redundant_widespread)

png("./figs/Growth_forms/Total_vs_redundant_species.png", width = 800, height = 500)
ggplot(data = new_df, 
       mapping = aes(x = Biome, fill = GROWTHFORM_STD, 
                     y = ifelse(test = Tmnt == "Total", 
                                yes = -prop, no = prop))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values=col_GF)+
  labs(y = "%") +
  coord_flip()+
  geom_hline(yintercept=0)
dev.off()


# Trait distribution: Redundant and widespread species --------------------

## Total distribution
cols=wes_palette("Chevalier")[c(1,3,4)]

dir.create("./figs/Trait_dist")

pdf("./figs/Trait_dist/Total_height_dis.pdf")
Traits_Biome_Di_Ri %>%
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Height)),
                      #calc_ecdf = TRUE,
                      rel_min_height = 0.01,
                      scale=0.9,na.rm = TRUE,alpha=0.8,fill= cols[1])+
  ggtitle("Total")+
  scale_y_discrete(labels=gsub("_"," ", levels(Traits_Biome_Di_Ri$Biome)))+
  xlim(-5,5)
dev.off()

pdf("./figs/Trait_dist/Red_Wides_height_dis.pdf", width = 5.5)
Traits_Biome_Di_Ri %>% 
  filter(Ri<=0.5 & DiScale < 0.2) %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = log(Height)),
                       #calc_ecdf = TRUE,
                       rel_min_height = 0.01,
                       scale=0.9,na.rm = TRUE,alpha=0.8,fill= cols[1])+
  ggtitle("Redundant & widespread species")+
  scale_y_discrete(labels=NULL)+
  xlim(-5,5)+
  ylab(" ")
dev.off()


pdf("./figs/Trait_dist/Total_height_dis_GF.pdf",width = 8)
Traits_Biome_Di_Ri %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges2(aes(x = log(Height), fill = paste(Biome, GROWTHFORM_GEN)),
                       #calc_ecdf = TRUE,
                       rel_min_height = 0.01,
                       stat = "density",
                       scale=1,na.rm = TRUE,alpha = .8, color = "white")+
  xlim(-3,4.5)+
  scale_fill_cyclical(values = cols)+
  theme_ridges() +
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_text(vjust = 0))+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "No information", "Woody"),
                      name = "Growth form", guide = "legend")+
  scale_y_discrete(labels=gsub("_"," ", levels(Traits_Biome_Di_Ri$Biome)))+
  ggtitle("Total")

dev.off()

pdf("./figs/Trait_dist/Red_Wides_height_dis_GF.pdf",width = 8)
Traits_Biome_Di_Ri %>% 
  filter(Ri<=0.5 & DiScale < 0.2) %>% 
  ggplot(aes(x=log(Height), y=Biome, height=..density..)) +
  geom_density_ridges2(aes(x = log(Height), fill = paste(Biome, GROWTHFORM_GEN)),
                       #calc_ecdf = TRUE,
                       rel_min_height = 0.01,
                       stat = "density",
                       scale=1,na.rm = TRUE,alpha = .8, color = "white")+
  xlim(-3,4.5)+
  scale_fill_cyclical(values = cols)+
  theme_ridges() +
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_text(vjust = 0))+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "No information", "Woody"),
                      name = "Growth form", guide = "legend")+
  ggtitle("Redundant & widespread species")+
  scale_y_discrete(labels=gsub("_"," ", levels(Traits_Biome_Di_Ri$Biome)))

dev.off()