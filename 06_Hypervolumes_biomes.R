# libraries ---------------------------------------------------------------
library(hypervolume)
library(raster)
library(rworldmap)
library(tidyverse)
library(visdat)
library(RColorBrewer)
library(foreach)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

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

biomes_names<-unique(cells_biomes$biomes)

Biomes_Abun_sp<-foreach(i=1:length(biomes_names),.combine=rbind)%do%
{
  
  indx<-which(row.names(biome_PAbs_matrix)==biomes_names[i])
  biome_abun<-colSums(biome_PAbs_matrix[indx,])
  
}
row.names(Biomes_Abun_sp)<-biomes_names

## Create dataframe with Species and biomes

Biomes_df<-Biomes_Abun_sp
Biomes_df[which(Biomes_df>0)]<-1

Biomes_df_clean<-
  Biomes_df%>%
  as.matrix %>% 
  t()%>%
  as.data.frame()%>%
  mutate(species=colnames(Biomes_df))%>%
  gather(key="Biome",value="Presence",-species) %>% 
  filter(Presence>0) %>% 
  dplyr::select(-Presence)


## Merging traits and distribution dataframe
Traits_phylo$species<-rownames(Traits_phylo)

Biome_traits_df<-merge(Traits_phylo, Biomes_df_clean)

# Calculate hypervolumes -------------------------------

# Scalling variables
Biome_traits_df$logseed_mass<-log(Biome_traits_df$seed_mass+1)


## Include Growth form
indx<-match(Biome_traits_df$species, Growth_form$SPECIES_STD)
Biome_traits_df$GROWTHFORM_STD<-Growth_form$GROWTHFORM_STD[indx]


woody<-c("Tree","Liana","Shrub","Woody epiphyte")
Biome_traits_df$GROWTHFORM_GEN<-ifelse(Biome_traits_df$GROWTHFORM_STD%in%woody,"woody","herbaceous")
Biome_traits_df$GROWTHFORM_GEN[which(is.na(Biome_traits_df$GROWTHFORM_STD))]<-NA

#Selecting and Scalling variables
scale_fun<-function(x){scale(x, center=TRUE, scale=TRUE)}


Hyper_input<-
  foreach(i=1:length(biome_name)) %do%
  {
    herbs<-Biome_traits_df %>% 
      filter(Biome==biome_name[i] & GROWTHFORM_GEN=="herbaceous") %>% 
      select(logseed_mass, wood_density, height, leaf_area) %>% 
      apply(2,scale_fun)
    
    woody<-Biome_traits_df %>% 
      filter(Biome==biome_name[i] & GROWTHFORM_GEN=="woody") %>% 
      select(logseed_mass, wood_density, height, leaf_area) %>% 
      apply(2,scale_fun)
    
    tmp<-list(woody=woody, herbs=herbs)
    
}

names(Hyper_input)<-biome_name

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


##plotting hypevolumes
#Log transform amphibian hypervolume
plot(Hyper_final_2)

plot(Hyper_final,show.3d=TRUE,plot.3d.axes.id=2:4,cex.random=3,cex.data=6,
     show.legend=TRUE,point.alpha.min=0.5,point.dark.factor=1)


plot(hypervolume_join(moist, savanna),num.points.max.random=6000,contour.lwd=1.5,colors=c(brewer.pal(n=2,"Set1"))
     ,show.legend=FALSE)
legend("bottomleft",legend = c("moist","Savanna"),text.col=c(brewer.pal(n=2,"Set1")),bty="n",cex=1.1,text.font=2)
