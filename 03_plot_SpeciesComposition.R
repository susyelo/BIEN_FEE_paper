# libraries ---------------------------------------------------------------
library(raster)
library(tidyverse)
library(foreach)
library(fuzzySim)
library(circlize)
library(RColorBrewer)
library(wesanderson)
library(dendextend)

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


# 3. Shapefiles
biome_shp<-shapefile("./data/processed/Olson_processed/Biomes_olson_projected.shp")

## Include cell number using the row and col numbers as reference
spPresence$cells<-cellFromRowCol(r_Total_Rich,spPresence$Y, spPresence$X)


## Folders
#dir.create("./figs/species_composition")

###############
# PROCEDURE
###############

# 1. Extract species for each biome --------------------------------------------
## Base raster
r_base<-r_Total_Rich
values(r_base)<-1:ncell(r_base)
names(r_base)<-"cell"

p <- rasterToPolygons(r_base)

Cells_biomes<-
  foreach(i=1:length(biome_shp), .combine = rbind)%do%{
    
    print(paste("Extract cells from",biome_shp$biomes[i]))
    biome_tmp<-biome_shp[i,]
    
    v <-as.data.frame(over(biome_tmp,p,returnList = TRUE))
    v$biome<-biome_tmp$biomes
    
    v
  }

spPresence_biome<-merge(spPresence, Cells_biomes, by.x="cells", by.y="cell")
saveRDS(spPresence_biome, file="./outputs/spPresence_biomes_all.rds")


## Number of cells per species in each biome
cells_in_sp<-spPresence_biome %>% 
  group_by(Species,biome) %>% 
  summarise(N_cells=n_distinct(cells)) %>% 
  group_by(Species) %>% 
  mutate(Total_cells=sum(N_cells), prop_cells=N_cells/sum(N_cells)) %>% 
  mutate(max_prop=max(prop_cells))

saveRDS(cells_in_sp, file="./outputs/spPresence_cell_prop_biomes_all.rds")


# 2. Species list for each biome ------------------------------------------

# 2.1 Total numbr of species
Total_sp_list<-tapply(cells_in_sp$Species,cells_in_sp$biome,unique)

# 2.2 Species with highest proportion of their ranges in each biome
Wides_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==max_prop)

Wides_sp_list<-tapply(Wides_sp$Species,Wides_sp$biome,unique)

# 2.3 Endemics for each biome
Endemics_sp<-cells_in_sp %>% 
  dplyr::filter(prop_cells==1)

Endemics_sp_list<-tapply(Endemics_sp$Species,Endemics_sp$biome,unique)

# 2.4 Proportion of endemics in each biome
total_n<-unlist(lapply(Total_sp_list,length))
endemics_n<-unlist(lapply(Endemics_sp_list,length))
prop_endemics<-round(endemics_n/total_n,3)*100


# 3. Create similarity matrix ---------------------------------------------
## Create a loop to calculate the similarity (number of species shared among biomes)
biome_richness<-Total_sp_list

spSimilarity<-foreach(i=1:length(biome_richness), .combine='cbind') %:%
  foreach(j=1:length(biome_richness), .combine='c') %do% {
    length(intersect(biome_richness[[i]],biome_richness[[j]]))
  }

colnames(spSimilarity)<-names(biome_richness)
rownames(spSimilarity)<-names(biome_richness)

## Double check that the numbers are correct
#diag(spSimilarity)==unlist(lapply(biome_richness, n_distinct))
#spSimilarity/diag(spSimilarity)

# 4. Chordplot of similarities --------------------------------------------

# 4.1 Species composition among biomes using all the species
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

pdf("./figs/species_composition/Total_similarity_biomes_withEndemics.pdf")
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


# 4.2 Species composition using Dominant species (those that occupy most of their ranges in each biome)
# This is to give some direction of shared species among biome 
spSimilarity_Wides<-foreach(i=1:length(Wides_sp_list), .combine='rbind') %:%
  foreach(j=1:length(Total_sp_list), .combine='rbind') %do% {
    
    N_sp=length(intersect(Wides_sp_list[[i]],Total_sp_list[[j]]))
    df<-data.frame(from=names(Wides_sp_list)[i], 
                   to=names(Total_sp_list)[j],
                   Sp_shared=N_sp)
    df
  }

spSimilarity_ma<-matrix(data = spSimilarity_Wides$Sp_shared, nrow = n_distinct(spSimilarity_Wides$from), ncol = n_distinct(spSimilarity_Wides$to), byrow = FALSE,
            dimnames = NULL)

rownames(spSimilarity_ma)<-unique(spSimilarity_Wides$from)
colnames(spSimilarity_ma)<-unique(spSimilarity_Wides$to)

## Print file to include into the supplementary information
write.csv(spSimilarity_ma,"./supp_info/Shared_species_matrix.csv")

diag(spSimilarity_ma)<-0
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))

# Order by richness
indx<-rev(order(Wides_sp_total))
spSimilarity_ma<-spSimilarity_ma[indx,indx]

colnames(spSimilarity_ma)<-c("Moist","Trop Dry",
                            "Xeric","Trop Grass",
                            "Savannas",
                            "Coniferous","Temp Mixed",
                            "Temp Grass","Mediterranean","Taiga","Tundra")


## Proportion of widespread species
total_n<-unlist(lapply(Total_sp_list,length))
Wides_sp_total<-unlist(lapply(Wides_sp_list,length))
prop_widespread<-round(Wides_sp_total/total_n,2)*100


colnames(spSimilarity_ma)<-paste(colnames(spSimilarity_ma),", ", prop_widespread[indx],"%", sep="")
rownames(spSimilarity_ma)<-colnames(spSimilarity_ma)

pdf("./figs/species_composition/Total_similarity_biomes_DominantSp.pdf",width = 8)
par(mar=c(0, 0, 0, 0))
chordDiagram(spSimilarity_ma,column.col = col,
             grid.col =col, directional = -1, 
             direction.type = c("diffHeight", "arrows"),link.largest.ontop=TRUE,
             annotationTrack = c("grid"),link.arr.length = 0.2,link.arr.type = "big.arrow",
             preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.8)
}, bg.border = NA)
dev.off()



# 5. Species dissimilarity among biomes --------------------------------------
Similarity_sp_biomes<-function(sp_list){
  
  biome_similarity <- foreach(i = 1:length(sp_list), .combine='rbind') %:%
    foreach(j = 1:length(Total_sp_list), .combine='rbind') %do% {
      
      sp_intersection = length(intersect(sp_list[[i]],sp_list[[j]]))
      biome1 <- length(sp_list[[i]])
      biome2 <- length(sp_list[[j]])
      
      sorensen <- 2 * sp_intersection/(biome1 + biome2)
      sorensen_df<-data.frame(sorensen=sorensen, 
                             from=names(sp_list)[i], 
                             to=names(sp_list)[j])
      sorensen_df
    }
  
  Similarity_ma<-matrix(data = biome_similarity$sorensen, 
                        nrow = n_distinct(biome_similarity$from), 
                        ncol = n_distinct(biome_similarity$to), byrow = FALSE,
                          dimnames = NULL)
  
  rownames(Similarity_ma)<-unique(biome_similarity$from)
  colnames(Similarity_ma)<-unique(biome_similarity$to)
  
  Similarity_ma
}
  
## Dissimilarity all species  
Total_similarity<-Similarity_sp_biomes(Total_sp_list)
fit_total_sim <-hclust(as.dist(1-Total_similarity))

labels(fit_total_sim)<-c("Xeric", "Trop grass", "Savannas", "Trop Dry", "Moist", "Taiga","Tundra","Mediterranean",
  "Coniferous","Temp Grass","Temp Mixed")

dend_total<-
  fit_total_sim %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[2]) %>% 
  set("branches_lwd", 4)  %>% 
  set("labels_cex", 1.5)

pdf("./figs/species_composition/species_composition_cluster_allsp.pdf", height = 9.4, width = 9.1)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()



## Dissimilarity with dominant species
Dominant_similarity<-Similarity_sp_biomes(Wides_sp_list)
fit_Dominant_similarity <-hclust(as.dist(1-Dominant_similarity))

## Check the order first
#labels(fit_Dominant_similarity)
labels(fit_Dominant_similarity)<-c("Taiga","Tundra", "Moist", "Trop grass", 
                                   "Savannas", "Trop Dry", "Xeric",
                                   "Coniferous","Mediterranean",
                                   "Temp Grass","Temp Mixed")
  
dend_dom<-
  fit_Dominant_similarity %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[3]) %>% 
  set("branches_lwd", 4) %>% 
  set("labels_cex", 1.5)

pdf("./figs/species_composition/species_composition_cluster_Dominant_sp.pdf", height = 10, width = 9.1)
circlize_dendrogram(dend_dom,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()
