# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(foreach)
library(fuzzySim)
library(circlize)
library(RColorBrewer)
library(wesanderson)

# functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

# data --------------------------------------------------------------------
# 1. Presence of species in cells
## I am still not sure the order of the variables in the table
# TODO: Check this! 
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv",header = FALSE,
                     col.names = c("Species","Y","X"))


#2. Total_richness raster
r_Total_Rich<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

## Include NA to the cells that have not information on them
r_Total_Rich[r_Total_Rich==0]<-NA

#spplot(r_Total_Rich)

# 3. Shapefiles
load("./data/base/Danilo_data/biomes_shp4susy.RData")

## Include cell number using the row and col numbers as reference
spPresence$cells<-cellFromRowCol(r_Total_Rich,spPresence$Y, spPresence$X)


# Extract species for each biome --------------------------------------------
# Ignoring montane TODO: do it with better polygons
biome_poly<-list(moist,tropical.mixed,
                 savanna,grasslands,
                 dry,xeric,mediterranean,
                 temperate.mixed,coniferous,prairies,taiga,tundra)
  
names(biome_poly)<-c("moist","tropical.mixed",
                     "savanna","grasslands",
                     "dry","xeric","mediterranean",
                     "temperate.mixed","coniferous","prairies","taiga","tundra")


## Include biome classification in each cell 
spPresence$biomes<-NA
for (i in 1:length(biome_poly)){
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,biome_poly[[i]]))
  spPresence$biomes[spPresence$cells%in%cells_tmp]<-names(biome_poly)[i]
  
}

save(spPresence, file="./outputs/spPresence_biomes_all.RData")

## Presence/absence matrix of biomes
tmp<-spPresence %>% 
  select(Species, biomes) %>%
  splist2presabs(sites.col = "biomes", sp.col = "Species")

save(tmp, file="./outputs/Biome_ALL_Sp_matrix.RData")

## square matrix of pair-wise similarities among biomes
# biome.sim.mat<-simMat(tmp[,-1], method = "Jaccard",upper=FALSE)



biome_richness<-foreach(i=1:length(biome_poly))  %do% {
  
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,biome_poly[[i]]))
  sp_list_tmp<-unique(spPresence$Species[spPresence$cells%in%cells_tmp])
  
}
names(biome_richness)<-names(biome_poly)


## Create similarity matrix
## Create a loop to calculate the similarity (number of species shared among biomes)

spSimilarity<-foreach(i=1:length(biome_richness), .combine='cbind') %:%
  foreach(j=1:length(biome_richness), .combine='c') %do% {
    length(intersect(biome_richness[[i]],biome_richness[[j]]))
  }

colnames(spSimilarity)<-names(biome_richness)
rownames(spSimilarity)<-names(biome_richness)

## Double check that the numbers are correct
#diag(spSimilarity)==unlist(lapply(biome_richness, n_distinct))


spSimilarity/diag(spSimilarity)
# Chordplot of similarities -----------------------------------------------

col=c(wes_palette("Darjeeling",6,type="continuous"),
      wes_palette("Cavalcanti",6,type="continuous"))

spSimilarity_1<-spSimilarity
diag(spSimilarity_1)<-0
chordDiagram(spSimilarity_1, transparency = 0.25,
             grid.col =col,link.sort = TRUE,symmetric = TRUE)


colnames(spSimilarity_1)<-c("moist","tropical.mixed",
                            "savanna","grasslands",
                            "dry","xeric","MED",
                            "TEMP","CONI","PRA","TA","TU")

rownames(spSimilarity_1)<-colnames(spSimilarity_1)

pdf("./figs/Total_similarity_biomes.pdf")
chordDiagram(spSimilarity_1, grid.col =col,symmetric = TRUE,link.sort = TRUE,
             column.col = col)
dev.off()




# Plot with species with traits -------------------------------------------
trait_BIEN<-read.csv("./data/processed/BIEN_trait_GrowthForm.csv")
trait_BIEN$scrubbed_species_binomial<-gsub(" ","_",trait_BIEN$scrubbed_species_binomial)

spPresence_sub<-spPresence[spPresence$Species%in%trait_BIEN$scrubbed_species_binomial,]
spPresence_sub<-na.omit(spPresence_sub)


# Species list
biome_richness_sub<-foreach(i=1:length(biome_poly))  %do% {
  
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,biome_poly[[i]]))
  sp_list_tmp<-unique(spPresence_sub$Species[spPresence_sub$cells%in%cells_tmp])
  
}
names(biome_richness_sub)<-names(biome_poly)


## Create similarity matrix
## Create a loop to calculate the similarity (number of species shared among biomes)
spSimilarity_sub<-foreach(i=1:length(biome_richness_sub), .combine='cbind') %:%
  foreach(j=1:length(biome_richness_sub), .combine='c') %do% {
    length(intersect(biome_richness_sub[[i]],biome_richness_sub[[j]]))
  }

colnames(spSimilarity_sub)<-names(biome_richness_sub)
rownames(spSimilarity_sub)<-names(biome_richness_sub)


col=c(wes_palette("Darjeeling",6,type="continuous"),
      wes_palette("Cavalcanti",6,type="continuous"))

spSimilarity_sub1<-spSimilarity_sub
diag(spSimilarity_sub1)<-0

colnames(spSimilarity_sub1)<-c("moist","tropical.mixed",
                            "savanna","grasslands",
                            "dry","xeric","MED",
                            "TEMP","CONI","PRA","TA","TU")

rownames(spSimilarity_sub1)<-colnames(spSimilarity_sub1)


pdf("./figs/OnlyTrait_similarity_biomes.pdf")
chordDiagram(spSimilarity_sub1, grid.col =col,symmetric = TRUE,link.sort = TRUE,
             column.col = col)
dev.off()



# Plot species with growth form -------------------------------------------
Growth_form<-read.table("data/base/GrowthForm_Final.txt", header = TRUE)
Growth_form$SPECIES_STD<-gsub(" ","_",Growth_form$SPECIES_STD)

spPresence_sub<-spPresence[spPresence$Species%in%Growth_form$SPECIES_STD,]
spPresence_sub<-na.omit(spPresence_sub)


# Species list
biome_richness_sub<-foreach(i=1:length(biome_poly))  %do% {
  
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,biome_poly[[i]]))
  sp_list_tmp<-unique(spPresence_sub$Species[spPresence_sub$cells%in%cells_tmp])
  
}
names(biome_richness_sub)<-names(biome_poly)


## Create similarity matrix
## Create a loop to calculate the similarity (number of species shared among biomes)
spSimilarity_sub<-foreach(i=1:length(biome_richness_sub), .combine='cbind') %:%
  foreach(j=1:length(biome_richness_sub), .combine='c') %do% {
    length(intersect(biome_richness_sub[[i]],biome_richness_sub[[j]]))
  }

colnames(spSimilarity_sub)<-names(biome_richness_sub)
rownames(spSimilarity_sub)<-names(biome_richness_sub)


col=c(wes_palette("Darjeeling",6,type="continuous"),
      wes_palette("Cavalcanti",6,type="continuous"))

spSimilarity_sub1<-spSimilarity_sub
diag(spSimilarity_sub1)<-0

colnames(spSimilarity_sub1)<-c("moist","tropical.mixed",
                               "savanna","grasslands",
                               "dry","xeric","MED",
                               "TEMP","CONI","PRA","TA","TU")

rownames(spSimilarity_sub1)<-colnames(spSimilarity_sub1)


pdf("./figs/OnlyGrowth_similarity_biomes.pdf")
chordDiagram(spSimilarity_sub1, grid.col =col,symmetric = TRUE,link.sort = TRUE,
             column.col = col)
dev.off()





