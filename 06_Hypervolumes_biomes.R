# libraries ---------------------------------------------------------------
library(hypervolume)
library(tidyverse)
library(RColorBrewer)
library(foreach)
library(BIEN)
library(factoextra)

# Functions ---------------------------------------------------------------
source("./functions/BIEN2.0_RangeMaps_functions.R")

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
  



## Calculate biome's hypervolumes

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


### Redundant and widespread species hypervolumes

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
plot(
  hypervolume_join(
    Total_hypervol$Moist_Forest, 
    Total_hypervol$Dry_Forest
  ),
  contour.lwd=1.5,
  colors=c(brewer.pal(n=3,"Set1")),
  show.legend=TRUE
)



## Calculate Hypervolume similarity using Sorense's index

choices<-choose(length(biome_names),2) #x choose 2 possible pairs
combs<-combn(length(biome_names),2) # create those pairs

similarity<-NULL
A<-NULL
B<-NULL

for(i in 1:choices) {
  nums<-combs[,i]
  
  sim<-
    hypervolume_overlap_statistics(
      hypervolume_set(
        Total_hypervol[[nums[1]]], 
        Total_hypervol[[nums[2]]],
        check.memory = FALSE
      )
    )
  
  similarity<-rbind(similarity,sim)
  
  A[i]<-nums[1]
  B[i]<-nums[2]
}

###
Similarity_df<-data.frame(Biome1=biome_names[A],Biome2=biome_names[B])
Similarity_df<-cbind(Similarity_df,similarity)

tmp<-matrix(data = NA, nrow = 11, ncol = 11, byrow = FALSE,
            dimnames = NULL)

rownames(tmp)<-c(as.character(Similarity_df$Biome1[1]),as.character(unique(Similarity_df$Biome2)))
colnames(tmp)<-c(as.character(unique(Similarity_df$Biome1)),as.character(Similarity_df$Biome2[10]))

diag(tmp)<-1

tmp[lower.tri(tmp, diag = FALSE)]<-Similarity_df$sorensen

fit <-hclust(as.dist(1-tmp))
plot(fit)

## With the total values

similarity_hypervol<-function(list_hyper){

  choices<-choose(length(biome_names),2) #x choose 2 possible pairs
  combs<-combn(length(biome_names),2) # create those pairs
  
  similarity_total<-NULL
  A<-NULL
  B<-NULL
  
  for(i in 1:choices) {
    nums<-combs[,i]
    
    sim<-
      hypervolume_overlap_statistics(
        hypervolume_set(
          biomes_hypervolumes_total[[nums[1]]], 
          biomes_hypervolumes_total[[nums[2]]],
          check.memory = FALSE
        )
      )
    
    similarity_total<-rbind(similarity_total,sim)
    
    A[i]<-nums[1]
    B[i]<-nums[2]
  }

    return(similarity_total)  
}

###

Create_dist_matrix<-function(Similarity_df, index="sorensen"){
  
  tmp<-matrix(data = NA, nrow = 11, ncol = 11, byrow = FALSE,
              dimnames = NULL)
  
  rownames(tmp)<-c(as.character(Similarity_df$Biome1[1]),as.character(unique(Similarity_df$Biome2)))
  colnames(tmp)<-c(as.character(unique(Similarity_df$Biome1)),as.character(Similarity_df$Biome2[10]))
  
  diag(tmp)<-1
  tmp[lower.tri(tmp, diag = FALSE)]<-Similarity_df[[index]]
  
  return(tmp)

}

## Biome cluster using total species
Sim_total<-data.frame(Biome1=biome_names[A],Biome2=biome_names[B])
Sim_total<-cbind(Sim_total,similarity_total)

dist_total<-Create_dist_matrix(Sim_total)

fit_total <-hclust(as.dist(1-dist_total))
labels(fit_total)<-c("Moist","Trop_grass","Dry","Savannas","Taiga",
               "Tundra","Xeric_wood","Mediterranean","Coniferous","Temp_grass","Temp_mixed")

dend_total<-
  fit_total %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[1]) %>% 
  set("branches_lwd", 4)

#dir.create("./figs/hypervolumes_clusters")

pdf("./figs/hypervolumes_clusters/Total_Sorensen.pdf",height = 9.4, width = 9.1)
circlize_dendrogram(dend_total,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()

## Biome cluster using Redundant and widespread species----

Sim_red<-data.frame(Biome1=biome_names[A],Biome2=biome_names[B])
Sim_red<-cbind(Sim_red,similarity)

dist_red<-Create_dist_matrix(Sim_red)

fit_red <-hclust(as.dist(1-dist_red))
labels(fit_red)<-c("Moist","Trop_grass","Dry","Savannas","Taiga",
                   "Tundra","Xeric_wood","Coniferous","Temp_Mixed",
                   "Temp_grass","Mediterranean")


dend_red<-
  fit_red %>% 
  as.dendrogram() %>% 
  color_branches(1,col=wes_palette("Cavalcanti")[3]) %>% 
  set("branches_lwd", 4)


#dir.create("./figs/hypervolumes_clusters")

pdf("./figs/hypervolumes_clusters/Redundant_Sorensen.pdf", height = 9.4, width = 9.1)
par(mar=c(0, 0, 0, 0))
circlize_dendrogram(dend_red,dend_track_height = 0.7,labels_track_height = 0.2)
dev.off()
