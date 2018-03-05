# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(foreach)

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

library(circlize)
library(RColorBrewer)
library(wesanderson)

col=c(wes_palette("Darjeeling",6,type="continuous"),
      wes_palette("Cavalcanti",6,type="continuous"))


adj = matrix(sample(c(1, 0), 26**2, replace = TRUE, prob = c(1, 9)),
             nrow = 26, dimnames = list(LETTERS, LETTERS))
adj = ifelse(adj == 1, runif(26**2), 0)
chordDiagram(adj, transparency = 0.4, grid.col = "midnightblue",
             col = colorRamp2(seq(0, 1, 0.2), brewer.pal(6, "Blues")))

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



