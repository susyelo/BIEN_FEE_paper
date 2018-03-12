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

rescaleRas <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x,na.rm=TRUE)
  if(is.null(x.max)) x.max = max(x,na.rm=TRUE)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}


# data --------------------------------------------------------------------
# 1. Trait data frame
DF_trait_biome<-read.csv("./outputs/Df_traits_biomes_Sp_withRanges.csv")

#2. Presence and absence matrix of species
x=load("./data/processed/PresAbs_matrix_TraitRange.RData")
PAbs_matrix<-get(x)

#3. Richness raster as reference
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

#4. Biomes shapefiles
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

# Check raster and shapefiles are in the same projection
#proj4string(r_Total_Rich)==proj4string(biome_shp)


# Subsetting data ---------------------------------------------------------
# Selecting continuous trait variables and species 
DF_trait_biome_sub<-
  DF_trait_biome%>%
  dplyr::select(scrubbed_species_binomial,SLA,Seed_mass,Height,Leaf_P,Leaf_N,Wood_density) %>% 
  distinct

rownames(DF_trait_biome_sub)<-DF_trait_biome_sub$scrubbed_species_binomial

# Selecting species in the preAbs matrix that have at least one trait
indx<-which(colnames(PAbs_matrix)%in%DF_trait_biome_sub$scrubbed_species_binomial)
PAbs_matrix<-PAbs_matrix[,indx]


# Extract cells per biome -------------------------------------------------
## Include biome classification in each cell 
biome_name<-biome_shp$biomes

cells_biomes<-data.frame(cells=rownames(PAbs_matrix))
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
biome_PAbs_matrix<-PAbs_matrix
indx<-match(rownames(PAbs_matrix),cells_biomes$cells)
rownames(biome_PAbs_matrix)<-cells_biomes$biomes[indx]

biome_PAbs_matrix<-biome_PAbs_matrix[!is.na(rownames(biome_PAbs_matrix)),]



# Compute distance matrix of trait between each pair of species  ----------
Dist_matrix<-compute_dist_matrix(DF_trait_biome_sub[,-1],metric="euclidean",center = TRUE,
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

Biomes_di = distinctiveness_mod(Biomes_pabs_sp, Dist_matrix)


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
            sd_di=sd(Di))

summary(Biomes_di_clean$Di[which(Biomes_di_clean$Biome=="Moist_Forest")])

### Merging growth form data with distinctiveness data
df1<-
  DF_trait_biome %>% 
  dplyr::select(scrubbed_species_binomial,GROWTHFORM_GEN,GROWTHFORM_STD) %>% 
  distinct

Biome_di_growthForm<-merge(df1,Biomes_di_clean,
                           by.x="scrubbed_species_binomial",
                           by.y="species")

biome_name<-unique(Biome_di_growthForm$Biome)

### Distribution of Distinctiveness per biome ----
## Overall
cols=wes_palette("Chevalier")[c(1,3,4)]

Biome_di_growthForm %>% 
  ggplot(aes(x=Di, y=Biome, height=..density..)) +
  geom_density_ridges()


Biome_di_growthForm %>% 
  group_by(Biome) %>% 
  dplyr::summarise(mean_di=mean(Di),
                   sd_di=sd(Di), 
                   max=max(Di))

## There are around 94 species with extreme values
Biome_di_growthForm %>% 
  filter(Di>5) %>%
  group_by(Biome) %>% 
  dplyr::summarise(N_sp=n_distinct(scrubbed_species_binomial))

Biome_di_growthForm<-
  Biome_di_growthForm %>% 
  filter(Di<5)

Biome_di_growthForm %>% 
  ggplot(aes(x=Di, y=Biome, height=..density..)) +
  geom_density_ridges()

## By Growthform
Biome_di_growthForm %>% 
  ggplot(aes(x=Di, y=Biome, height=..density..)) +
  geom_density_ridges(aes(x = Di, fill = paste(Biome, GROWTHFORM_GEN)),
                      scale=2,na.rm = TRUE,alpha = .8, color = "white")+
  scale_fill_cyclical(values = cols,
                      labels = c("Herbaceous", "No information", "Woody"),
                      name = "Growth form", guide = "legend")+
  theme_ridges(grid = FALSE)



dist<-
foreach(i=1:length(biome_name), .combine = rbind) %do%{
  a<-Biome_di_growthForm %>% 
    filter(Biome==biome_name[i]) %>% 
    filter(Di>=quantile(Di, 0.7)) %>% 
    group_by(GROWTHFORM_STD) %>% 
    dplyr::summarise(N_sp=length(GROWTHFORM_STD)) %>% 
    mutate(Dist="Dist",Biome=biome_name[i],prop=round(N_sp/sum(N_sp)*100,2))
}

redunt<-
  foreach(i=1:length(biome_name), .combine = rbind) %do%{
    a<-Biome_di_growthForm %>% 
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

pdf("./figs/Growth_total.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = total_dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Total")
dev.off()



pdf("./figs/Growth_distinctive.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = dist,
                    stat="identity")+
  coord_flip() +
  ggtitle("Distinctive species")
dev.off()


pdf("./figs/Growth_redundant.pdf", width = 12)
ggplot() + geom_bar(aes(y = prop, x = Biome, fill = GROWTHFORM_STD), data = redunt,
                          stat="identity")+
  coord_flip() +
  ggtitle("Redundant species")
dev.off()  


ggplot(redunt) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(Biome), y=prop, fill=GROWTHFORM_STD),stat="identity") +
  coord_flip() 



