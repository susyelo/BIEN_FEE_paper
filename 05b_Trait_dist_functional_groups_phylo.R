# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)
library(sp)
library(foreach)
library(raster)
library(viridis)
library(ggridges)
library(wesanderson)


# Functions ---------------------------------------------------------------
source("./functions/funrar_functions_mod.R")
source("./functions/check_functions.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")

rescaleRas <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x,na.rm=TRUE)
  if(is.null(x.max)) x.max = max(x,na.rm=TRUE)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}


# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.table("./data/base/Danilo_data/traits4susy.csv", header = TRUE)

#2. Presence matrix of species
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv")
names(spPresence) = c("Species","Y","X")

#3. Richness raster as reference
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

#4. Biomes shapefiles
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

#5. Growth form
Growth_form<-read.table("./data/base/GrowthForm_Final.txt",header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)

# Check raster and shapefiles are in the same projection
#proj4string(r_Total_Rich)==proj4string(biome_shp)

# Presence absence matrix -------------------------------------------------
TraitSpecies <- unique(rownames(Traits_phylo))
spMatrix <- splistToMatrix(spPresence,TraitSpecies)

# Extract cells per biome
## Include biome classification in each cell 
biome_name<-biome_shp$biomes

cells_biomes<-data.frame(cells=rownames(spMatrix))
cells_biomes$biomes<-NA

for (i in 1:length(biome_name)){
  print(paste("Extracting",biome_name[i]))
  shp_tmp<-biome_shp[which(biome_shp$biomes==biome_name[i]),]
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,shp_tmp))
  cells_tmp<-paste("Cell",cells_tmp,sep="_")
  cells_biomes$biomes[cells_biomes$cells%in%cells_tmp]<-biome_name[i]
}

cells_biomes<-na.omit(cells_biomes)

## Biome presence and absence matrix
biome_PAbs_matrix<-spMatrix
indx<-match(rownames(spMatrix),cells_biomes$cells)
rownames(biome_PAbs_matrix)<-cells_biomes$biomes[indx]

# Filter cells without information
biome_PAbs_matrix<-biome_PAbs_matrix[!is.na(rownames(biome_PAbs_matrix)),]


# Compute distance matrix of trait between each pair of species  ----------
Dist_matrix<-compute_dist_matrix(Traits_phylo,metric="euclidean",center = TRUE,
                                 scale = TRUE) ## This can take a while



# Compute functional distinctiveness per biome ----------------------------

# Matrix with number of times the species is present in a grid per biome
# This could be use as relative abundance of species in each of the biomes
biomes_names<-unique(row.names(biome_PAbs_matrix))

Biomes_Abun_sp<-foreach(i=1:length(biomes_names),.combine=rbind)%do%
{
  
  indx<-which(row.names(biome_PAbs_matrix)==biomes_names[i])
  biome_abun<-colSums(biome_PAbs_matrix[indx,])
  
}
row.names(Biomes_Abun_sp)<-biomes_names

# Calculating relative abundance
Biome_relAbun<-make_relative(Biomes_Abun_sp)

# Calculating distinctiveness ----
Biomes_pabs_sp<-Biomes_Abun_sp
Biomes_pabs_sp[which(Biomes_pabs_sp>0)]<-1

Biomes_di = distinctiveness(Biomes_pabs_sp, Dist_matrix)

# Calculating distinctiveness using relative abundance (number of grids) of species ----
#Biomes_di_abundance = distinctiveness_mod(Biome_relAbun, Dist_matrix)

Biomes_di_clean<-
  Biomes_di%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_di))%>%
  gather(key="Biome",value="Di",-species) %>% 
  filter(!is.na(Di))


### Exploring data
Biomes_di_clean %>% 
  group_by(Biome) %>% 
  dplyr::summarise(mean_di=mean(Di),
            sd_di=sd(Di), max=max(Di))


### Include Growth form
indx<-match(Biomes_di_clean$species, Growth_form$SPECIES_STD)
Biomes_di_clean$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]

biome_name<-unique(Biomes_di_clean$Biome)

### Distribution of Distinctiveness per biome ----
## Overall

Biomes_di_clean %>% 
  ggplot(aes(x=Di, y=Biome, height=..density..)) +
  geom_density_ridges()


Biomes_di_clean %>% 
  group_by(Biome) %>% 
  dplyr::summarise(mean_di=mean(Di),
                   sd_di=sd(Di), 
                   max=max(Di))

## There are around 122 species with extreme values
Biomes_di_clean %>% 
  filter(Di>8) %>%
  group_by(Biome) %>% 
  dplyr::summarise(N_sp=n_distinct(species))

Biomes_di_clean<-
  Biomes_di_clean %>% 
  filter(Di<8)

Biomes_di_clean %>% 
  ggplot(aes(x=Di, y=Biome, height=..density..)) +
  geom_density_ridges()


dist<-
foreach(i=1:length(biome_name), .combine = rbind) %do%{
  a<-Biomes_di_clean %>% 
    filter(Biome==biome_name[i]) %>% 
    filter(Di>=quantile(Di, 0.7)) %>% 
    group_by(GROWTHFORM_STD) %>% 
    dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
    mutate(Dist="Dist",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,2))
}

redunt<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Biomes_di_clean %>% 
      filter(Biome==biome_name[i]) %>% 
      filter(Di<quantile(Di, 0.2)) %>% 
      group_by(GROWTHFORM_STD) %>% 
      dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
      mutate(Dist="Redun",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,2))
  }

dist<-
dist %>% 
  filter(GROWTHFORM_STD!="Aquatic")

redunt<-
  redunt %>% 
  filter(GROWTHFORM_STD!="Aquatic")



total_dist<-rbind(dist,redunt)

## Total growth form proportion

pdf("./figs/Growth_total_phylo.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = total_dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Total")
dev.off()



pdf("./figs/Growth_distinctive_phylo.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Distinctive species")
dev.off()


pdf("./figs/Growth_redundant_phylo.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = redunt,
                          stat="identity")+
  coord_flip() +
  ggtitle("Redundant species")
dev.off()  


ggplot(redunt) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(Biome), y=prop, fill=GROWTHFORM_STD),stat="identity") +
  coord_flip() 

