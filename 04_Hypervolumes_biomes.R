# libraries ---------------------------------------------------------------
library(hypervolumes)
library(raster)
library(rworldmap)

# data --------------------------------------------------------------------
# Olson shapefiles and raster files
load("./data/Danilo_data/biomes_shp4susy.RData")
load("./data/Olson_raster.RData")

occ_data<-readRDS("./data/2018_02_08_BIEN_traitOccurData.RData")

occ_data_geo<-occ_data
coordinates(occ_data_geo)<-cbind(occ_data_geo$longitude , occ_data_geo$latitude)

biomes_names<-levels(r_biomes)[[1]]

### Changing projections

proj4string(occ_data_geo)<-proj4string(r_biomes)

crs(occ_data_geo)<-crs(r_biomes)
occ_data_geo_trnsfrmd = spTransform(occ_data_geo,crs(r_biomes))


occ_data_geo$biomes<-extract(r_biomes,occ_data_geo,factors=TRUE)

match(occ_data_geo$biomes, biomes_names$ID)

plot(wrld_simpl,xlim=c(-179,-40),ylim=c(-20,80))
points(occ_data_geo_trnsfrmd[is.na(occ_data_geo_trnsfrmd$biomes),],col="red",pch=20)
