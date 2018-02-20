# libraries ---------------------------------------------------------------
library(letsR)


# data --------------------------------------------------------------------
#trait data
trait_df<-read.csv("./outputs/BIEN_trait_GrowthForm.csv")

# Occurrence data
occ_df<-readRDS("./data/2018_02_08_BIEN_traitOccurData.RData")

# Filter Occurrence data --------------------------------------------------
Occ_BIEN_sub<-
  occ_df %>% 
  filter(scrubbed_species_binomial %in% unique(trait_df$scrubbed_species_binomial))


# Create species richness map
xy<-data.frame(x=Occ_BIEN_sub$longitude, y=Occ_BIEN_sub$latitude)

# Presence/abscence raster of Solanum ---
BIEN_trait_grids<-lets.presab.points(xy, Occ_BIEN_sub$scrubbed_species_binomial, 
                                     resol=1,count = TRUE,
                                     show.matrix=TRUE)

BIEN_trait_grids$Richness_Raster[BIEN_trait_grids$Richness_Raster==0]<-NA


## Create a new raster
resol<-1
r<-raster(xmn = -167, xmx = -29,ymn = -57, ymx = 83,crs = CRS("+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
res(r) <- resol

## Find cells with attributes
cells <- cellFromXY(r, Occ_BIEN_sub[,c("longitude","latitude")])

tmp<-rasterize(xy,r,field=1)
species_data<-Occ_BIEN_sub[,c("longitude","latitude")]

pa.raster <- presence.absence.raster(mask.raster=r, species.data=species_data, raster.label="tmp")
plot(pa.raster, main=names(pa.raster))


## Ploting the richness of Solanum
e <- extent(-170, 180, -60, 75)
ras<-spplot(crop(log(Solanum_grids$Richness_Raster),e),
            sp.layout=spl_tmp,col.regions=colors, colorkey=list(labels=list(cex=2)))


