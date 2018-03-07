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
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

## Include cell number using the row and col numbers as reference
spPresence$cells<-cellFromRowCol(r_Total_Rich,spPresence$Y, spPresence$X)


# Extract species for each biome --------------------------------------------

biome_name<-biome_shp$biomes

## Include biome classification in each cell 
spPresence$biomes<-NA
for (i in 1:length(biome_name)){
  print(paste("Extracting",biome_name[i]))
  shp_tmp<-biome_shp[which(biome_shp$biomes==biome_name[i]),]
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,shp_tmp))
  spPresence$biomes[spPresence$cells%in%cells_tmp]<-biome_name[i]
}

save(spPresence, file="./outputs/spPresence_biomes_all.RData")

## Presence/absence matrix of biomes
#tmp<-spPresence %>% 
#  select(Species, biomes) %>%
#  splist2presabs(sites.col = "biomes", sp.col = "Species")

save(tmp, file="./outputs/Biome_ALL_Sp_matrix.RData")

## square matrix of pair-wise similarities among biomes
# biome.sim.mat<-simMat(tmp[,-1], method = "Jaccard",upper=FALSE)

biome_richness<-foreach(i=1:length(biome_name))  %do% {
  
  print(paste("Extracting",biome_name[i]))
  shp_tmp<-biome_shp[which(biome_shp$biomes==biome_name[i]),]
  cells_tmp<-unlist(cellFromPolygon(r_Total_Rich,shp_tmp))
  sp_list_tmp<-unique(spPresence$Species[spPresence$cells%in%cells_tmp])
  
}
names(biome_richness)<-biome_name

## Calculate the number of endemic species

spEndemics<-foreach(i=1:length(biome_richness), .combine='c') %do%{
  
  compare<-which(names(biome_richness)!=names(biome_richness)[i])
  N_endemics<-setdiff(biome_richness[[i]],unlist(biome_richness[compare]))
  n_distinct(N_endemics)
}
names(spEndemics)<-names(biome_richness) 

# Total number of species per biome
total_n<-unlist(lapply(biome_richness,n_distinct))

prop_endemics<-round(spEndemics/total_n,2)*100


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
#spSimilarity/diag(spSimilarity)
# Chordplot of similarities -----------------------------------------------

## order by richness 
indx<-rev(order(diag(spSimilarity)))
spSimilarity_1<-spSimilarity
spSimilarity_1<-spSimilarity_1[indx,indx]


col=c(wes_palette("Darjeeling",6,type="continuous"),
      wes_palette("Cavalcanti",5,type="continuous"))

diag(spSimilarity_1)<-0
colnames(spSimilarity_1)<-c("Moist","Dry",
                            "Xeric","Savannas",
                            "Trop Grass","Coniferous","Temp Mixed",
                            "Temp Grass","Mediterranean","Taiga","Tundra")


colnames(spSimilarity_1)<-paste(colnames(spSimilarity_1),", ", prop_endemics[indx],"%", sep="")

rownames(spSimilarity_1)<-colnames(spSimilarity_1)

pdf("./figs/Total_similarity_biomes.pdf")
chordDiagram(spSimilarity_1, grid.col =col,symmetric = TRUE,
             column.col = col)
dev.off()



# now, the image with rotated labels

pdf("./figs/Total_similarity_biomes_withEndemics.pdf")
par(mar=c(0, 0, 0, 0))
chordDiagram(spSimilarity_1, annotationTrack = "grid", preAllocateTracks = 1, grid.col =col,symmetric = TRUE,
             column.col = col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .2, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.7)
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()




