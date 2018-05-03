## Calculate hypervolumes for each Biome function
Biomes_hypervolume<-function(biome_dataframe, biome_names){
  
  require(foreach)
  hyper_list<-foreach(i=seq_along(biome_names))%do%{
    
    biome_df<-
      biome_dataframe %>% 
      dplyr::filter(Biome==biome_names[i])
    
    biome_hb<-hypervolume_gaussian(biome_df[,-1],name = as.character(biome_names[i]))
    
    biome_hb
  }
  
  names(hyper_list)<-biome_names
  hyper_list
}

## Sorensen similarity of a list of hypervolumes
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


