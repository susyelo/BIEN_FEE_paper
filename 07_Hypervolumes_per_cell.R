# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(foreach)
library(raster)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")
source("./functions/Biomes_hypervolumes_fun.R")

# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB_lambda.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

# 5. Presence matrix of species
spPresence<-read.csv("./data/base/BIEN_2_Ranges/presence100km.csv")
names(spPresence) = c("Species","Y","X")
TraitSpecies <- unique(Traits_phylo$species)
spMatrix_sub <- splistToMatrix(spPresence,TraitSpecies)

## Remove cells without species
indx<-which(rowSums(spMatrix_sub)!=0)
spMatrix_sub<-spMatrix_sub[indx,]

# Data manipulation -------------------------------------------------------
# 1. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)


## 3. Rename and order biomes
Traits_Biome_Di_Ri$Biome<-recode(Traits_Biome_Di_Ri$Biome,Moist_Forest="Moist",
                                 Savannas="Savannas",
                                 Tropical_Grasslands="Trop_Grass",
                                 Dry_Forest="Dry",
                                 Xeric_Woodlands="Xeric",
                                 Mediterranean_Woodlands="Mediterranean",
                                 Temperate_Grasslands="Temp_Grass",
                                 Temperate_Mixed="Temp_Mixed",
                                 Coniferous_Forests="Coniferous",
                                 Taiga="Taiga",
                                 Tundra="Tundra")

Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome,
                                 levels=c("Moist","Savannas","Trop_Grass",
                                          "Dry","Xeric","Mediterranean",
                                          "Temp_Grass","Temp_Mixed","Coniferous",
                                          "Taiga","Tundra"))

# Calculate hypervolumes -------------------------------

# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density)
Traits_Biome_Di_Ri$sqrtSLA<-sqrt(Traits_Biome_Di_Ri$SLA)

#Selecting and Scalling variables
Traits_Biome_Di_Ri<-
  Traits_Biome_Di_Ri %>%
  mutate(Scaled_logSeed_mass=as.numeric(scale(logseed_mass)),
         Scaled_logHeight=as.numeric((logHeight)),
         Scaled_SLA=as.numeric(scale(sqrtSLA)),
         Scaled_logWood_density=as.numeric(scale(logWoodDensity)),
         Scaled_Leaf_N=as.numeric(scale(Leaf_N)),
         Scaled_Leaf_P=as.numeric(scale(Leaf_P))
  )

Trait_df<-
  Traits_Biome_Di_Ri %>%
  dplyr::select(species,contains("Scaled"))


## Hypervolumes for each biomes cells using hypervolume box algorithm
cell_biomes_df<-readRDS("./outputs/spPresence_biomes_all.rds")

cell_biomes<-cell_biomes_df %>%
  filter(Species%in%Trait_df$species)

cell_biomes <-tapply(cell_biomes$cells, cell_biomes$biomes, unique)

cells_names<-as.character(as.vector(unlist(cell_biomes)))

cells_names <- paste("Cell", cells_names, sep="_")



hypervolume_per_cell <- function(cells_list, pb_matrix, trait_df){
  
  count <- 0
  
  foreach(i = cells_list, .combine = rbind)%do%{
    
    print(i)
    count <- count + 1
    x<-pb_matrix[i,]
    
    print(paste("Processing",count, "out of ",length(cells_list)))
    
    sp_names<-names(x[x > 0 & !is.na(x)])
    
    if (length(sp_names)>1){
      
      res<- tryCatch({
        cell_hyper<-trait_df %>%
          filter(species%in%sp_names) %>%
          dplyr::select(contains("Scaled")) %>%
          hypervolume_box(kde.bandwidth=0.5)
        
        cell_hyper@Volume
        
      },
      error = function(cond){
        message("Species with the same trait values")
        return(NA)
      })
      
    }else{
      res = NA
    }
    
    tmp_df <- data.frame(cell = i, vol = res)
    
    tmp_df
  }
  
}

## Running the function 
system.time(
hyper_cells_box <- hypervolume_per_cell(cells_list = cells_names, 
                            pb_matrix = spMatrix_sub, 
                            trait_df = Trait_df)
)

write_rds(hyper_cells_box, "./outputs/07_hypervol_cell_SDM_BW0.5.rds")

## 
write_rds(hyper_cells_box, "./outputs/07_hypervol_cell_SDM_501_end_box.rds")

hyper_cells_box_1<-readRDS("outputs/07_hypervol_cell_SDM_501_end_box.rds")
hyper_cells_box_2<-readRDS("outputs/07_hypervol_cell_SDM_1:500_box.rds")

hyper_cells_box<- rbind(hyper_cells_box_2,hyper_cells_box_1)

hyper_cells_box$cell<-as.numeric(gsub("Cell_","",hyper_cells_box$cell))

indx<-match(hyper_cells_box$cell,cell_biomes_df$cells)
hyper_cells_box$biomes<-cell_biomes_df$biomes[indx]

hyper_cells_box$biomes<-factor(hyper_cells_box$biomes,
                             levels=c("Moist_Forest","Savannas","Tropical_Grasslands",
                                      "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                      "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                      "Taiga","Tundra"))

hyper_cells_box$biomes<-recode(hyper_cells_box$biomes,Moist_Forest="Moist",
                             Savannas="Savannas",
                             Tropical_Grasslands="Trop_Grass",
                             Dry_Forest="Dry",
                             Xeric_Woodlands="Xeric",
                             Mediterranean_Woodlands="Mediterranean",
                             Temperate_Grasslands="Temp_Grass",
                             Temperate_Mixed="Temp_Mixed",
                             Coniferous_Forests="Coniferous",
                             Taiga="Taiga",
                             Tundra="Tundra")

library(wesanderson)

pdf("./figs/07_Hypervolume_cells_SDM_box.pdf", width=10)
ggplot(data=hyper_cells_box,aes(x=biomes,y=vol)) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("SD"^"6")))
dev.off()


# Add richness info

spMatrix_sub_copy <-spMatrix_sub
spMatrix_sub_copy[which(spMatrix_sub_copy>0)]<-1

cell_richness <- rowSums(spMatrix_sub_copy)

names(cell_richness) <- gsub("Cell_", "", names(cell_richness))

indx<-match(hyper_cells_box$cell,names(cell_richness))
hyper_cells_box$Richness<-cell_richness[indx]


library(ggpmisc)
library(ggpubr)

biomes_to_plot <- c("Dry","Xeric","Mediterranean")

hyper_cells_box$logRich <- log(hyper_cells_box$Richness)

tmp_df <- hyper_cells_box %>% 
  filter(biomes%in%biomes_to_plot & logRich > 4)


png("./figs/07_Hypervolume_cells_SDM_box_DMX.png", width=700)
ggscatterhist(data = tmp_df, x = "logRich", y = "vol",
              color = "biomes", size = 3, alpha = 0.6,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              margin.plot = "boxplot",
              ggtheme = theme_bw())
dev.off()


## Moist, Temperate, and Taiga
biomes_to_plot <- c("Moist","Temp_Mixed","Taiga")

tmp_df <- hyper_cells_box %>% 
  filter(biomes%in%biomes_to_plot & logRich > 4)


png("./figs/07_Hypervolume_cells_SDM_box_MoisTempTaig.png", width=700)
ggscatterhist(data = tmp_df, x = "logRich", y = "vol",
              color = "biomes", size = 3, alpha = 0.6,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              margin.plot = "boxplot",
              ggtheme = theme_bw())
dev.off()




### 




### Taking only the 10% of the cells per each biomes
cell_biomes_df<-readRDS("./outputs/spPresence_biomes_all.rds")

cell_biomes<-cell_biomes_df %>%
  filter(Species%in%Trait_df$species)

cell_biomes <-tapply(cell_biomes$cells, cell_biomes$biomes, unique)

### Taking only the 20% of the cells per each biomes
Random_cells<-
  lapply(cell_biomes,
         function(x)
           sample(x,length(x)*.20)
  )

cells_names<-as.character(as.vector(unlist(Random_cells)))

cells_names <- paste("Cell", cells_names, sep="_")


Tmp<-NULL
count <- 0

system.time(

for (i in cells_names)
{
  print(i)
  count <- count + 1
  x<-spMatrix_sub[i,]

  print(paste("Processing",count, "out of ",length(cells_names)))

  sp_names<-names(x[x > 0 & !is.na(x)])

  #if(length(sp_names)>100){

   # sample_sp<-sample(sp_names,100)
  #}else{
   # sample_sp<-sp_names
  #}

  if (length(sp_names)>1){

    res<- tryCatch({
      cell_hyper<-Trait_df %>%
        filter(species%in%sp_names) %>%
        dplyr::select(contains("Scaled")) %>%
        hypervolume_gaussian()

      cell_hyper@Volume

    },
    error = function(cond){
      message("Species with the same trait values")
      return(NA)
    })

  }else{

    res=NA

  }

  tmp_df <- data.frame(cell = i, vol = res)

  Tmp<-rbind(Tmp,tmp_df)

  write_rds(Tmp,"./outputs/07_Hypervolume_sp_sample_box.rds")
}
)

cell_hyper_df <- Tmp
cell_hyper_df$cell<-as.numeric(gsub("Cell_","",cell_hyper_df$cell))

indx<-match(cell_hyper_df$cell,cell_biomes_df$cells)
cell_hyper_df$biomes<-cell_biomes_df$biomes[indx]

cell_hyper_df$biomes<-factor(cell_hyper_df$biomes,
                             levels=c("Moist_Forest","Savannas","Tropical_Grasslands",
                                      "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                      "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                      "Taiga","Tundra"))

cell_hyper_df$biomes<-recode(cell_hyper_df$biomes,Moist_Forest="Moist",
                                 Savannas="Savannas",
                                 Tropical_Grasslands="Trop_Grass",
                                 Dry_Forest="Dry",
                                 Xeric_Woodlands="Xeric",
                                 Mediterranean_Woodlands="Mediterranean",
                                 Temperate_Grasslands="Temp_Grass",
                                 Temperate_Mixed="Temp_Mixed",
                                 Coniferous_Forests="Coniferous",
                                 Taiga="Taiga",
                                 Tundra="Tundra")

library(wesanderson)

pdf("./figs/07_Hypervolume_cells_spNewHypervolume_box.pdf", width=10)
ggplot(data=cell_hyper_df,aes(x=biomes,y=vol)) +
  geom_boxplot()+
  geom_jitter(alpha=0.5,color=wes_palette("Cavalcanti1")[4])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab(expression(paste("SD"^"6")))
dev.off()


# Add richness info
spMatrix_sub[which(spMatrix_sub>0)]<-1

cell_richness <- rowSums(spMatrix_sub)

names(cell_richness) <- gsub("Cell_", "", names(cell_richness))

indx<-match(cell_hyper_df$cell,names(cell_richness))
cell_hyper_df$Richness<-cell_richness[indx]



library(ggpmisc)
library(ggpubr)

biomes_to_plot <- c("Dry","Xeric","Mediterranean")

cell_hyper_df$logRich <- log(cell_hyper_df$Richness)

tmp_df <- cell_hyper_df %>% 
  filter(biomes%in%biomes_to_plot & logRich > 4)

ggscatterhist(data = tmp_df, x = "logRich", y = "vol",
              color = "biomes", size = 3, alpha = 0.6,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              margin.plot = "boxplot",
              ggtheme = theme_bw())


Total_richness <- raster("./data/base/BIEN_2_Ranges/richness100km.tif")

tmp_df$TotallogRich <- log(Total_richness[tmp_df$cell])
## Loop of 50 times


### Taking only the 10% of the cells per each biomes
cell_biomes_df<-readRDS("./outputs/spPresence_biomes_all.rds")

cell_biomes<-cell_biomes_df %>%
  filter(Species%in%Trait_df$species)

cell_biomes <-tapply(cell_biomes$cells, cell_biomes$biomes, unique)


for (j in 1:50){

  Random_cells<-
    lapply(cell_biomes,
           function(x)
             sample(x,length(x)*.10)
    )

  Random_cells<-lapply(Random_cells, function(x)paste("Cell",x,sep="_"))

  cells_names<-as.vector(unlist(Random_cells))

  Tmp<-NULL

  for (i in cells_names)
  {
    x<-spMatrix_sub[i,]

    print(paste("Processing",length(Tmp)))

    cell_names<-names(x[x > 0 & !is.na(x)])

    if(length(cell_names)>10){

      sample_sp<-sample(cell_names,10)
    }else{
      sample_sp<-cell_names
    }

    if (length(cell_names)>1){

      cell_hyper<-Trait_df %>%
        filter(species%in%sample_sp) %>%
        dplyr::select(contains("Scaled")) %>%
        hypervolume_gaussian()

      res<-cell_hyper@Volume

    }else{

      res=0

    }
    Tmp<-c(Tmp,res)
    names(Tmp)<-cell_names
    write_rds(Tmp,"outputs/Hypervolume_sp_sample_gaussian50_Itera.rds")
  }

}
