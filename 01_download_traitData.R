# LIBS --------------------------------------------------------------------
library(BIEN)
library(tidyverse)
library(maptools)
library(rworldmap)
library(cleangeo)


# Download data -----------------------------------------------------------
spp <- BIEN_list_all()$species                              # Species list
traits <- BIEN_trait_species(spp)                           # Trait list

occur_spp<-BIEN_occurrence_species(spp, 
                                   natives.only=TRUE, 
                                   only.new.world=TRUE)     # Occurrence data only New world

# Filter data -----------------------------------------------------------
# Remove info without coordinates
occur_spp_Geo<-
  occur_spp%>%
  filter(!is.na(latitude) & !is.na(longitude))


# Write downloaded data ---------------------------------------------------
name_file<-paste(format(Sys.time(), "%Y_%m_%d"), "BIEN_trait", sep="_")

# Traits
write.csv(traits,
          file=paste("./data/",name_file,"data.csv", sep=""))

# Occurrences
saveRDS(occur_spp_Geo, file=paste("./data/",name_file,"OccurData.RData", sep=""))
