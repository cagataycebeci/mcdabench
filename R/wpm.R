# Weighted Product Method (WPM)
wpm <- function(dmatrix, bcvec, weights, normethod="vector", tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  
   nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = normethod) 
   
   wmatrix <- sweep(nmatrix, 2, weights, "^")
  
   p <- apply(wmatrix, 1, prod)
 
   ranking <- rank(-p, ties.method=tiesmethod)
   names(ranking) <- rownames(dmatrix)

   results <- list(
     nmatrix = nmatrix, 
     weighted_matrix = wmatrix, 
     p = p, 
     rank=ranking
   )
   
   return(results)
}
