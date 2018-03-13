# library -----------------------------------------------------------------
library(raster)
library(rgdal)
library(rgeos)
library(tidyverse)
library(tmap)
library(wesanderson)

# data --------------------------------------------------------------------
# 1. Biome polygons
download.file("http://maps.tnc.org/files/shp/terr-ecoregions-TNC.zip",
              destfile = "./data/base/wwf_olson/olson_biomes.zip")

### The unzip function does not work properly in windows so I need to set a temporary working directory
setwd("./data/base/wwf_olson/")
unzip("olson_biomes.zip")
setwd("../../../")


biomes_shp <- shapefile("./data/base/wwf_olson/tnc_terr_ecoregions.shp")


# Choosing New world polygons
inx<-which(biomes_shp$WWF_REALM2=="Nearctic" | biomes_shp$WWF_REALM2=="Neotropic")
biome_NW<-biomes_shp[inx,]


#2. Raster with the referenced projection (total richness raster)
total_richness<-raster("./data/base/BIEN_2_Ranges/richness100km.tif")

# Dissolving Ecoregion polygons -------------------------------------------
# Ensure shapefile row.names and polygon IDs are sensible
row.names(biome_NW) <- row.names(biome_NW@data)
biome_NW <- spChFIDs(biome_NW, row.names(biome_NW))

# Create new variable with the names we are using in BIEN 
biome_NW$biomes<-biome_NW$WWF_MHTNAM
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Temperate Broadleaf and Mixed Forests")]<-"Temperate_Mixed"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Tropical and Subtropical Coniferous Forests")]<-"Dry_Forest"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Temperate Grasslands, Savannas and Shrublands")]<-"Temperate_Grasslands"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Tropical and Subtropical Dry Broadleaf Forests")]<-"Dry_Forest"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Deserts and Xeric Shrublands")]<-"Xeric_Woodlands"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Tropical and Subtropical Moist Broadleaf Forests")]<-"Moist_Forest"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Flooded Grasslands and Savannas")]<-"Savannas"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Tropical and Subtropical Grasslands, Savannas and Shrublands")]<-"Savannas"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Montane Grasslands and Shrublands")]<-"Tropical_Grasslands"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Mediterranean Forests, Woodlands and Scrub")]<-"Mediterranean_Woodlands"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Temperate Conifer Forests")]<-"Coniferous_Forests"
biome_NW$biomes[which(biome_NW$WWF_MHTNAM=="Boreal Forests/Taiga")]<-"Taiga"

# Extract the data you want (the larger geography)
biome_poly <- gUnaryUnion(biome_NW, id = biome_NW@data$biomes)

# And add the data back in
lu<-row.names(biome_poly)
lu <- as.data.frame(lu)
colnames(lu) <- "biomes" 
row.names(lu)<-lu$biomes
 
biome_poly <- SpatialPolygonsDataFrame(biome_poly, lu)

# Filter some biomes ------------------------------------------------------
exclude_biomes<-c("Inland Water","Rock and Ice","Mangroves")
biome_poly_sub<-biome_poly[which(biome_poly$biomes%in%exclude_biomes==FALSE),]
biome_poly_sub$biomes<-droplevels(biome_poly_sub$biomes)

# Change shapefile projections --------------------------------------------
crs_ref<-crs(total_richness)
biomes_shp_proj <- spTransform(biome_poly_sub, CRSobj = crs_ref)

biomes_shp_proj$biomes<-factor(biomes_shp_proj$biomes, 
                               levels=c("Moist_Forest","Savannas","Tropical_Grasslands",
                                        "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                        "Temperate_Grasslands", "Coniferous_Forests",
                                        "Temperate_Mixed",
                                        "Taiga","Tundra"))




# Plot biomes shapefile ---------------------------------------------------
# manually set the colors for the plot!
bio_colors <- c('#33a02c','#b2df8a', '#ffff99','#fdbf6f','#ff7f00','#e31a1c','#fb9a99','#1f78b4','#a6cee3',
                '#cab2d6','#6a3d9a')

# plot using new colors
cols_wes=c(wes_palette("Darjeeling",6,type="continuous"),
           wes_palette("Cavalcanti",5,type="continuous"))


names_tmp<-gsub("_"," ",levels(biomes_shp_proj$biomes))

### Dataframe of colours for other analysis
Colours_biomes<-data.frame(biome=levels(biomes_shp_proj$biomes), color=bio_colors)
write.csv(Colours_biomes, "./outputs/colours_biomes.csv")

pdf("./figs/Biomes_shpfile1.pdf")
spplot(biomes_shp_proj,col.regions=bio_colors,
       pretty=T, names.attr=names_tmp)
dev.off()

pdf("./figs/Biomes_shpfile2.pdf")
qtm(biomes_shp_proj, fill="biomes", fill.style="fixed",
    fill.labels=names_tmp,
    fill.palette=bio_colors, 
    fill.title="Biomes")
dev.off()

pdf("./figs/Biomes_shpfile3.pdf")
qtm(biomes_shp_proj, fill="biomes", fill.style="fixed",
    fill.labels=names_tmp,
    fill.palette=cols_wes)
dev.off()


# Write shapefile ---------------------------------------------------------
writeOGR(obj=biomes_shp_proj, dsn="./data/processed/Olson_processed", layer="Biomes_olson_projected", driver="ESRI Shapefile")

