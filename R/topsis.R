topsis <- function(dmatrix, bcvec, weights, normethod = "maxmin", tiesmethod="average") {
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
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
   
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = normethod) 
 
  wmatrix <- sweep(nmatrix, MARGIN = 2, STATS = weights, FUN = "*")
   
  pis <- apply(wmatrix, 2, max)
  nis <- apply(wmatrix, 2, min)
  Dp <- apply(wmatrix, 1, function(x) {sqrt(sum((x - pis)^2))})
  Dm <- apply(wmatrix, 1, function(x) {sqrt(sum((x - nis)^2))})
  p <- Dm / (Dm + Dp)
  
  ranking <- rank(-p, ties.method=tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(
    normalized_matrix = nmatrix, 
    weighted_matrix = wmatrix, 
    nis = nis, 
    pis = pis, 
    Dm = Dm, 
    Dp = Dp, 
    p = p, 
    rank = ranking
 )
    
  return(results)
}

