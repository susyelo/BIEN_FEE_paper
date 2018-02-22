############
# spFromRowCol
# This function takes row and column numbers and Jim's big presence matrix
# and returns a character vector containing the species that occur in the
# given cell

# Arguments:
#  d = Jim's Presence Matrix
#  row = the row number to query
#  col = the column number to query

# Return: A character vector containing species names

spFromRowCol = function(d, row, col)
{
  index = which(d$Y == row & d$X == col)
  return(d$Species[index])
}
############

############
# spFromCell
# This function takes a cell index and Jim's big presence matrix and
# returns a character vector containing the species that occur in the
# given cell
# Cell index is counted like this:
# 1  2  3  4
# 5  6  7  8
# 9 10 11 12

# Arguments:
#  d = Jim's Presence Matrix
#  cell = cell number to query
#  r = (Optional) the raster to provide the geographic system. If blank, assume the raster is Jim's 100km raster

# Return: A character vector containing species names

spFromCell = function(d, cell, r=NULL)
{
  #If the raster is not given, assume it is the same as Jim's richness raster
  if(is.null(r))
  {
    r = raster(nrows = 146, ncols = 103, xmn = -5261554, xmx = 5038446, ymn = -7434988, ymx = 7165012, crs = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }
  
  row = rowFromCell(r,cell)
  col = colFromCell(r,cell)
  index = which(d$Y == row & d$X == col)
  return(d$Species[index])
}
############

############
# spListToMatrix
# This function takes the big matrix of presence data and a character vector
# containing species names. It returns a matrix with one column per given
# species name and one row per cell in the raster, with rows ordered in cell
# order (see above). Values in the matrix are 0 or 1 for absence or presence,
# or NA if the given species name did not match any names in Jim's matrix

# Arguments:
#  d = Jim's Presence Matrix
#  sp = A character vector containing species names
#  r = (Optional) the raster to provide the geographic system. If blank, assume the raster is Jim's 100km raster

# Return: A presence/absence matrix

splistToMatrix = function(d,sp,r = NULL)
{
  #If the raster is not given, assume it is the same as Jim's richness raster
  if(is.null(r))
  {
    r = raster(nrows = 146, ncols = 103, xmn = -5261554, xmx = 5038446, ymn = -7434988, ymx = 7165012, crs = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }
  
  #out will be the presence/absence matrix
  out = matrix(0,nrow = ncell(r),ncol = length(sp))
  d2 = d[which(d$Species %in% sp),]
  
  rowind = cellFromRowCol(r,d2$Y,d2$X)
  colind = match(d2$Species,sp)
  out[cbind(rowind,colind)] = 1
  
  out[,which(apply(out,2,sum)==0)] = NA
  rownames(out) = paste("Cell",1:ncell(r),sep ="_")
  colnames(out) = sp
  return(out)
}
############

############
# splistToRichness
# This function takes the presence matrix and a vector of species names
# and creates a raster of richness counts. This is very similar to spListToMatrix,
# but is easier on RAM, because #it only keeps track of a vector of species
# richness, not the whole presence/absence matrix.

# If a species name is given that does not match a name in the data, a warning
# is given.

# It returns the results as a raster, and optionally writes the
# raster to a file, if filename is given. If the file already exists, it will
# not be overwritten, and an error will be returned.

# Arguments:
#  d = Jim's Presence Matrix
#  sp = A character vector containing species names
#  r = (Optional) the raster to provide the geographic system. If blank, assume the raster is Jim's 100km raster

# Return: A richness raster (also, can write a file as a side effect)

splistToRichness = function(d,sp,r = NULL,filename = NULL)
{
  #If the raster is not given, assume it is the same as Jim's richness raster
  if(is.null(r))
  {
    r = raster(nrows = 146, ncols = 103, xmn = -5261554, xmx = 5038446, ymn = -7434988, ymx = 7165012, crs = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }
  d2 = d[which(d$Species %in% sp),]
  
  #create a richness vector
  out = rep(0,ncell(r))
  rowind = cellFromRowCol(r,d2$Y,d2$X)
  ta = table(rowind)
  
  out[as.numeric(names(ta))] = ta
  
  noMatch = sp[which(!sp%in%d2$Species)]
  #If there were unmatched species, print a warning
  if(length(noMatch) > 0)
  {
    print("Warning - there were unmatched species names. They were:")
    print(noMatch)
  }
  
  #Wrap the vector of richness values into a raster
  outRaster = setValues(r,out)
  
  #If a file name was given, write outRaster to that file (without overwriting)
  if(!is.null(filename))
  {
    writeRaster(outRaster,filename)
  }
  return(outRaster)
}
############

############
#genusToSplist
# This function takes the species name information from Jim's presence
# matrix, and given a character genus name, returns a character vector
# of species in that genus

#Note: I would like to extend this to also enable extraction by families,
# but do not currently have the taxonomic information to do that

# Arguments:
#  d = Jim's Presence Matrix
#  genus = The genus name

# Return: A character vector of species names

genusToSplist = function(d,genus)
{
  usp = unique(d$Species)
  splits = strsplit(usp,"_")
  splitgen = unlist(lapply(splits,function(i) {i[[1]][1]}))
  index = which(splitgen == genus)
  return(usp[index])
}

