# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)
library(factoextra)
library(dendextend)
library(wesanderson)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

similarity_hypervol<-function(list_hyper){
  
  choices<-choose(length(names(list_hyper)),2) #x choose 2 possible pairs
  combs<-combn(length(names(list_hyper)),2) # create those pairs
  
  similarity_total<-NULL
  A<-NULL
  B<-NULL
  
  for(i in 1:choices) {
    nums<-combs[,i]
    
    sim<-
      hypervolume_overlap_statistics(
        hypervolume_set(
          list_hyper[[nums[1]]], 
          list_hyper[[nums[2]]],
          check.memory = FALSE
        )
      )
    
    similarity_total<-rbind(similarity_total,sim)
    
    A[i]<-nums[1]
    B[i]<-nums[2]
    
    
  }
  
  Similarity_df<-data.frame(Biome1=names(list_hyper)[A],Biome2=names(list_hyper)[B])
  Similarity_df<-cbind(Similarity_df,similarity_total)
  
  tmp<-matrix(data = NA, nrow = n_distinct(list_hyper), ncol = n_distinct(list_hyper), byrow = FALSE,
              dimnames = NULL)
  
  rownames(tmp)<-names(list_hyper)
  colnames(tmp)<-names(list_hyper)
  
  diag(tmp)<-1
  
  tmp[lower.tri(tmp, diag = FALSE)]<-Similarity_df$sorensen
  
  return(tmp)  
}


# data --------------------------------------------------------------------
# 1. Trait data frame
Traits_phylo<-read.csv("./data/processed/traits_ALLMB.csv")

# 2. Values of distinctiveness and Restrictedness for species per biome
Biome_Di_Ri<-read.csv("./outputs/Biome_Di_Ri_phylo.csv", row.names = 1)

# 3. Merging data frames
Traits_Biome_Di_Ri<-merge(Biome_Di_Ri,Traits_phylo)

Traits_Biome_Di_Ri$Biome<-factor(Traits_Biome_Di_Ri$Biome, 
                                 levels=rev(c("Moist_Forest","Savannas","Tropical_Grasslands",
                                              "Dry_Forest","Xeric_Woodlands","Mediterranean_Woodlands",
                                              "Temperate_Grasslands","Temperate_Mixed","Coniferous_Forests",
                                              "Taiga","Tundra")))

# Calculate hypervolumes -------------------------------

# Transforming and Scaling variables
Traits_Biome_Di_Ri$logseed_mass<-log(Traits_Biome_Di_Ri$Seed_mass+1)
Traits_Biome_Di_Ri$logHeight<-log(Traits_Biome_Di_Ri$Height+1)
Traits_Biome_Di_Ri$logWoodDensity<-log(Traits_Biome_Di_Ri$Wood_density+1)

#Selecting and Scalling variables
Traits_Biome_Di_Ri<-
Traits_Biome_Di_Ri %>% 
  group_by(Biome) %>% 
  mutate(Scaled_logWood_density=scale(logWoodDensity),
         Scaled_Leaf_N=scale(Leaf_N),
         Scaled_SLA=scale(SLA),
         Scaled_Leaf_P=scale(Leaf_P),
         Scaled_logSeed_mass=scale(logseed_mass),
         Scaled_logHeight=scale(logHeight))
  



## Calculate hypervolumes for each Biome function
Biomes_hypervolume<-function(biome_dataframe, biome_names){
  
  hyper_list<-foreach(i=seq_along(biome_names))%do%{
    
    biome_df<-
      biome_dataframe %>% 
      dplyr::filter(Biome==biome_names[i]) %>% 
      dplyr::select(contains("Scaled"))
    
    biome_hb<-hypervolume_gaussian(biome_df[,-1],name = as.character(biome_names[i]))
    
    biome_hb
  }
  
  names(hyper_list)<-biome_names
  hyper_list
}


# Redundant and widespread species hypervolumes ---------------------------

## Ordering traits
Traits_Biome_Di_Ri<-
  Traits_Biome_Di_Ri %>% 
  select(species:logWoodDensity,
         Scaled_logHeight,Scaled_logWood_density,Scaled_logSeed_mass,
         Scaled_SLA,Scaled_Leaf_N,Scaled_Leaf_P)

biome_names<-unique(Traits_Biome_Di_Ri$Biome)

biome_ReduntWides<-
  Traits_Biome_Di_Ri %>% 
  dplyr::filter(Ri<=0.5 & DiScale < 0.2)

Redun_Wides_hypervol<-Biomes_hypervolume(biome_ReduntWides,biome_names)
saveRDS(Redun_Wides_hypervol, "./outputs/ReduntWides_hypervolumes.rds")

## Total hypervolumes
Total_hypervol<-Biomes_hypervolume(Traits_Biome_Di_Ri,biome_names)
saveRDS(Total_hypervol, "./outputs/Total_hypervolumes.rds")

## Plot hypervolumes per category
png("./figs/hypervolumes_clusters/Total_Moist_Dry_Savanna.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Total_hypervol$Moist_Forest, 
    Total_hypervol$Dry_Forest,
    Total_hypervol$Savannas
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/Total_Moist_Temperated_Mixed_Conifers.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Total_hypervol$Moist_Forest, 
    Total_hypervol$Temperate_Mixed,
    Total_hypervol$Coniferous_Forests
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/RedWides_Moist_Dry_Savanna.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist_Forest, 
    Redun_Wides_hypervol$Dry_Forest,
    Redun_Wides_hypervol$Savannas
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()

png("./figs/hypervolumes_clusters/RedWides_Moist_Temperated_Mixed_Conifers.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Redun_Wides_hypervol$Moist_Forest, 
    Redun_Wides_hypervol$Temperate_Mixed,
    Redun_Wides_hypervol$Coniferous_Forests
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)
dev.off()


## Calculate Hypervolume similarity using Sorense's index -----
## With Redundant species
redun_Sim<-similarity_hypervol(Redun_Wides_hypervol)
fit_red <-hclust(as.dist(1-redun_Sim))
labels(fit_red)<-c("Moist","Trop_grass","Dry","Savannas","Taiga",
                   "Tundra","Xeric_wood","Coniferous","Temp_Mixed",
                   "Temp_grass","Mediterranean")

dend_red<-
  fit_red %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[3]) %>% 
  set("branches_lwd", 4)

## With the total species
Total_Sim<-similarity_hypervol(Total_hypervol)
fit_total <-hclust(as.dist(1-Total_Sim))
labels(fit_total)<-c("Moist","Trop_grass","Dry","Savannas","Taiga",
               "Tundra","Xeric_wood","Mediterranean","Coniferous","Temp_grass","Temp_mixed")

dend_total<-
  fit_total %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[1]) %>% 
  set("branches_lwd", 4)

#dir.create("./figs/hypervolumes_clusters")
pdf("./figs/hypervolumes_clusters/Redundant_Sorensen.pdf", height = 9.4, width = 9.1)
circlize_dendrogram(dend_red,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

pdf("./figs/hypervolumes_clusters/Total_Sorensen.pdf",height = 9.4, width = 9.1)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

# Hypervolumes for climatic categories ------------------------------------
tropical<-c("Moist_Forest","Dry_Forest","Dry_Forest","Tropical_Grasslands","Savannas")
temperate<-c("Temperate_Grasslands","Coniferous_Forests","Temperate_Mixed","Mediterranean_Woodlands")
cold<-c("Taiga","Tundra")

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%tropical,"Tropical",NA)
Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%temperate,"Temperate",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes<-ifelse(Traits_Biome_Di_Ri$Biome%in%cold,"Cold",
                                      Traits_Biome_Di_Ri$ClimBiomes)

Traits_Biome_Di_Ri$ClimBiomes[which(Traits_Biome_Di_Ri$Biome=="Xeric_Woodlands")]<-"Xeric"


## The function needs a variable called "Biome" 
## so I created a temporal dataframe with the Biome variables using the climate classification
Traits_Biome_Di_Ri_tmp<-Traits_Biome_Di_Ri
Traits_Biome_Di_Ri_tmp$Biome<-Traits_Biome_Di_Ri_tmp$ClimBiomes

biome_names<-unique(Traits_Biome_Di_Ri_tmp$Biome)
Climatic_hypervol<-Biomes_hypervolume(Traits_Biome_Di_Ri_tmp,biome_names)
saveRDS(Climatic_hypervol,"./outputs/Climatic_hypervolumes.rds")


Climatic_Sim<-similarity_hypervol(Climatic_hypervol)
fit_Climatic <-hclust(as.dist(1-Climatic_Sim))


png("./figs/hypervolumes_clusters/Climatic_hypervolumes_total.png",width = 600, height = 600)
plot(
  hypervolume_join(
    Climatic_hypervol$Tropical, 
    Climatic_hypervol$Xeric,
    Climatic_hypervol$Temperate,  
    Climatic_hypervol$Cold
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=4,"Set1")),
  show.legend=TRUE,
  cex.data=1,cex.axis=1,cex.names=1,cex.legend=2,
  show.random=FALSE
)
dev.off()
