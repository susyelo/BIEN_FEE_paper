# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)
library(sp)
library(foreach)
library(raster)
library(viridis)
library(ggridges)
library(wesanderson)
library(lattice)


# Functions ---------------------------------------------------------------
source("./functions/check_functions.R")
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/trait_distribution_functions.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

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

# Presence absence matrix of species in each biome -------------------------------------------------
TraitSpecies <- unique(Traits_phylo$species)
spMatrix_sub <- splistToMatrix(spPresence,TraitSpecies)

Biomes_pabs<-spMatrix_biomes(spMatrix = spMatrix_sub,
                             biome_shp = biome_shp,
                             raster_ref = r_Total_Rich)

# Compute distance matrix of trait between each pair of species  ----------
rownames(Traits_phylo)<-Traits_phylo$species

# Select traits to calculate the distances among species
traits<-c("Wood_density","Leaf_N","SLA","Seed_mass","Height","Leaf_P")


traits<-c("Wood_density","Seed_mass")

Dist_matrix<-compute_dist_matrix(Traits_phylo[,traits],metric="euclidean",
                                 center = TRUE,
                                 scale = TRUE) ## This can take a while



# Compute functional distinctiveness per biome ----------------------------

# Calculating relative abundance
Biome_relAbun<-make_relative(Biomes_pabs$Biomes_Abun)

Biomes_di = distinctiveness(Biome_relAbun, Dist_matrix)

Biomes_di_clean<-
  Biomes_di%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_di))%>%
  gather(key="Biome",value="Di",-species) %>% 
  filter(!is.na(Di))


biome_names=biome_shp$biomes

# Compute functional restrictiness per biome ----------------------------
## A 0 value indicates that the focal species is present in all the sites.
rest_species<-
  foreach(i=1:length(biome_names), .combine = rbind)%do%
  {
    indx<-which(rownames(Biomes_pabs$Biomes_pabs_cells)==biome_name[i])
    biome_PAbs_tmp<-Biomes_pabs$Biomes_pabs_cells[indx,]
    rest_species<-restrictedness(pres_matrix = biome_PAbs_tmp)
    rest_species$Biome<-biome_name[i]
    
    rest_species
  }

## Exclude all that are not present in any biome cell (which are not in Biome_Di)
# Usually species restricted to just one cell which is not classified in any of the biomes
rest_species<-
  rest_species %>% 
  filter(species%in%unique(Biomes_di_clean$species))


## Scaling Di values per biome
Biomes_di_clean<-
  Biomes_di_clean %>% 
  group_by(Biome) %>% 
  mutate(DiScale=rescaleRas01(Di))

## Merge Di and Ri
Biome_Di_Ri<-merge(Biomes_di_clean, rest_species)

Biome_Di_Ri$FunDi<-(Biome_Di_Ri$DiScale+Biome_Di_Ri$Ri)/2
Biome_Di_Ri$FunDi2<-Biome_Di_Ri$Di*Biome_Di_Ri$Ri

write.csv(Biome_Di_Ri, "./outputs/Biome_Di_Ri_phylo.csv")


## Heatmaps
# Total headmap

Biome_Di_Ri$bin_Di<-cut(Biome_Di_Ri$DiScale, breaks = 10)
Biome_Di_Ri$bin_Ri<-cut(Biome_Di_Ri$Ri, breaks = 10)

get_matrix<-function(df=Biome_Di_Ri, Biome_name="Total")
{
  if(Biome_name=="Total"){
    tmp<-df
  }else{
    tmp<-df %>% 
      dplyr::filter(Biome==Biome_name) %>% 
      droplevels()
  }
  data<-table(tmp$bin_Ri,tmp$bin_Di)
}

biome_name=unique(Biome_Di_Ri$Biome)

for(i in 1:length(biome_name)){
  png(paste("./figs/Di_Ri_heatmaps/Heatmap_", biome_name[i],".png",sep=""))
  levelplot(log(get_matrix(df=Biome_Di_Ri,Biome_name=biome_name[i])+1),
            col.regions = heat.colors(100)[length(heat.colors(100)):1], 
            xlab="Ri",ylab="Di",ylim=c(levels(Biome_Di_Ri$bin_Di)),
            xlim=c(levels(Biome_Di_Ri$bin_Ri)),
            scales=list(x=list(rot=90)), main=biome_name[i])
  dev.off()
}



png("./figs/Taiga_Ri_Di.png")
levelplot(log(data+1),
          col.regions = heat.colors(100)[length(heat.colors(100)):1], 
          xlab="Ri",ylab="Di")
dev.off()


### Exploring data
Biome_Di_Ri %>% 
  group_by(Biome) %>% 
  dplyr::summarise(median_di=median(FunDi),
            sd_di=sd(FunDi), max=max(FunDi))


### Include Growth form
indx<-match(Biome_Di_Ri$species, Growth_form$SPECIES_STD)
Biome_Di_Ri$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]

biome_name<-unique(Biome_Di_Ri$Biome)

### Distribution of Distinctiveness per biome ----
## Overall


Biome_Di_Ri %>% 
  ggplot(aes(x=FunDi, y=Biome, height=..density..)) +
  geom_density_ridges()


## There are around 122 species with extreme values
dist<-
foreach(i=1:length(biome_name), .combine = rbind) %do%{
  a<-Biome_Di_Ri %>% 
    filter(Biome==biome_name[i]) %>% 
    filter(FunDi>=quantile(FunDi, 0.7)) %>% 
    group_by(GROWTHFORM_STD) %>% 
    dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
    mutate(Dist="Dist",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,2))
}

redunt<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Biome_Di_Ri %>% 
      filter(Biome==biome_name[i]) %>% 
      filter(FunDi<quantile(FunDi, 0.2)) %>% 
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

pdf("./figs/Growth_total_phyloSeed.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = total_dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Total")
dev.off()



pdf("./figs/Growth_distinctive_phyloSeed.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Distinctive species")
dev.off()


pdf("./figs/Growth_redundant_phyloSeed.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = redunt,
                          stat="identity")+
  coord_flip() +
  ggtitle("Redundant species")
dev.off()  


ggplot(redunt) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(Biome), y=prop, fill=GROWTHFORM_STD),stat="identity") +
  coord_flip() 


