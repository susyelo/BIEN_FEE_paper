library(Matrix)
require(funrar)

is_relative = function(given_obj, abund = NULL) {
  
  is_rel = FALSE
  
  if (is.matrix(given_obj) | is(given_obj, "sparseMatrix")) {
    values = na.omit(unique(as.vector(given_obj)))
  } else if (is.data.frame(given_obj) & !is.null(abund) & is.character(abund)) {
    values = na.omit(unique(given_obj[[abund]]))
  }
  
  max_val = max(values)
  min_val = min(values)
  
  if (max_val <= 1 & min_val >= 0) {
    is_rel = TRUE
  }
  
  return(is_rel)
}

# Modify function to being able to multiply the presence absence matrix with distance matrix that have NAs
# Nas appear when two species are compared but do not have enough trait information to calculate the distance with other species
distinctiveness_mod<-function (pres_matrix, dist_matrix) 
{
  full_matrix_checks(pres_matrix, dist_matrix)
  common = species_in_common(pres_matrix, dist_matrix)
  pres_matrix = pres_matrix[, common, drop = FALSE]
  dist_matrix = dist_matrix[common, common]
  if (!is_relative(pres_matrix)) {
    warning("Provided object may not contain relative abundances nor ", 
            "presence-absence\n", "Have a look at the make_relative() function if it is the case")
  }
  
  mat1m<- as(pres_matrix,"dgCMatrix")
  dist_matrix[is.na(dist_matrix)] <- 0
  mat2m<- as(dist_matrix,"dgCMatrix")
  
  index_matrix=tcrossprod(mat1m,mat2m)
  
  if (requireNamespace("Matrix", quietly = TRUE) & is(pres_matrix, 
                                                      "sparseMatrix")) {
    index_matrix[Matrix::which(pres_matrix == 0)] = NA
    total_sites = Matrix::rowSums(pres_matrix)
  }
  else {
    index_matrix[which(pres_matrix == 0)] = NA
    total_sites = rowSums(pres_matrix)
  }
  denom_matrix = apply(pres_matrix, 2, function(x) total_sites - 
                         x)
  index_matrix = index_matrix/denom_matrix
  if (any(vapply(index_matrix, function(x) is.nan(x), logical(1)))) {
    warning("Some communities had a single species in them\n", 
            "Computed value assigned to 'NaN'")
  }
  dimnames(index_matrix) = dimnames(pres_matrix)
  return(index_matrix)
}