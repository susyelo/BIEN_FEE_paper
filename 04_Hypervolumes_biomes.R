
# libraries ---------------------------------------------------------------
library(hypervolumes)


# data --------------------------------------------------------------------
load("./data/Danilo_data/biomes_shp4susy.RData")
biomes_grid<-read.table("./data/Danilo_data/biomes.csv")

load("./outputs/BIEN_traits_PresAbs1degree.RData")

occ_data<-readRDS("./data/2018_02_08_BIEN_traitOccurData.RData")

occ_data_geo<-occ_data


coordinates(occ_data_geo)<-cbind(occ_data_geo$longitude , occ_data_geo$latitude)


tmp<-extract(occ_data_geo,dry@data)
