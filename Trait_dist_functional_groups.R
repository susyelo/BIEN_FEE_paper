# library -----------------------------------------------------------------
library(funrar)
library(tidyverse)
library(sp)
library(foreach)
library(raster)


# Functions ---------------------------------------------------------------
source("./functions/funrar_functions_mod.R")
source("./functions/check_functions.R")

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

# Calculating distinctiveness using relative abundance (number of grids) of species ----
Biomes_di = distinctiveness_mod(Biome_relAbun, Dist_matrix)
Biomes_di_tmp = distinctiveness(Biome_relAbun, Dist_matrix)

tmp<-colSums(Biomes_di,na.rm = TRUE)
length(which(tmp>0))
## Check this file for further analysis
# /Users/echeverrialondono1/Git_Repos/BIEN_project/Functional_Div/4__Dominants_Subordinates_Dist.R
