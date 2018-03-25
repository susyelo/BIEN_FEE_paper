

rescaleRas01 <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x,na.rm=TRUE)
  if(is.null(x.max)) x.max = max(x,na.rm=TRUE)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}


fromCell_to_Biomes<-function(biome_shp, raster_ref)
{
  require(foreach)
  # Extract the names of the biomes from the shapefile
  biomes_names= biome_shp$biomes
  
  foreach(i=1:length(biomes_names), .combine = rbind)%do%{
    
    print(paste("Extracting",biomes_names[i]))
    shp_tmp<-biome_shp[which(biome_shp$biomes==biomes_names[i]),]
    cells_tmp<-unlist(cellFromPolygon(raster_ref,shp_tmp))
    cells_tmp<-paste("Cell",cells_tmp,sep="_")
    df_cell_biomes<-data.frame(cells=cells_tmp,biomes=biomes_names[i])
    
  }
  
}


spMatrix_biomes<-function(spMatrix, biome_shp, raster_ref){
  
  # spMatrix= Species presence absence matrix
  # biomes shapefile
  # Reference raster for projections and position of cells
  cells_biomes<-fromCell_to_Biomes(biome_shp = biome_shp, 
                                   raster_ref = raster_ref)
  
  
  ## Species presence and absence in biomes matrix
  biome_PAbs_matrix<-spMatrix
  indx<-match(rownames(spMatrix),cells_biomes$cells)
  rownames(biome_PAbs_matrix)<-cells_biomes$biomes[indx]
  biome_PAbs_matrix<-biome_PAbs_matrix[!is.na(rownames(biome_PAbs_matrix)),] # Remove cells without species
  
  
  # 1. Matrix with number of times the species is present in a grid per biome
  biomes_names<-unique(row.names(biome_PAbs_matrix))

  Biomes_Abun_sp<-foreach(i=1:length(biomes_names),.combine=rbind)%do%
  {
    
    indx<-which(row.names(biome_PAbs_matrix)==biomes_names[i])
    biome_abun<-colSums(biome_PAbs_matrix[indx,])
    
  }
  row.names(Biomes_Abun_sp)<-biomes_names
  
  #2. Final presence/absence matrix species vs biomes
  Biomes_pabs_sp<-Biomes_Abun_sp
  Biomes_pabs_sp[which(Biomes_pabs_sp>0)]<-1
  
  
  Biomes_pabs_final<-list(Biomes_Abun=Biomes_Abun_sp, 
                          Biomes_pabs=Biomes_pabs_sp, 
                          Biomes_pabs_cells=biome_PAbs_matrix)
}




## Di vs RI heatmaps
Di_Ri_heatmaps<-function(Biome_Di_Ri, xvar, yvar, breaks=10, Biome_toPlot, 
                         col.regions= heat.colors(100)[length(heat.colors(100)):1]){
  
  
  Biome_Di_Ri$bin_Di<-cut(yvar, breaks = breaks,dig.lab = 1,include.lowest = TRUE)
  Biome_Di_Ri$bin_Ri<-cut(xvar, breaks = breaks,dig.lab = 1,include.lowest = TRUE)
  
  Biome_Di_Ri_tmp<-
    Biome_Di_Ri %>% 
    dplyr::filter(Biome==Biome_toPlot) %>% 
    droplevels()
  
  data<-table(Biome_Di_Ri_tmp$bin_Ri,Biome_Di_Ri_tmp$bin_Di)
  
  levelplot(log(data),
            col.regions = col.regions, 
            xlab="Ri",ylab="Di",
            ylim=as.character(seq(0.05,1,by=0.1)),
            xlim=as.character(seq(0.05,1,by=0.1)),
            scales=list(x=list(rot=90)), main=Biome_toPlot)
  
  
}

