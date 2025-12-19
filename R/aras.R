# ARAS: Additive Ratio Assessment
aras <- function(dmatrix, bcvec, weights, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("Benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and 
       have the same length as the number of criteria.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
 
  if (!is.numeric(weights) || length(weights) != m) {
    stop("Weights must be a numeric vector with the same length as the number of criteria.")
  } 

  # Normalized Decision Matrix
  nmatrix <- calcnormal(dmatrix, bcvec, type="ratio")

  # Weighted Normalized Decision Matrix
  wmatrix <- sweep(nmatrix, 2, weights, "*")

  # Optimality function value (Si)
  Si <- rowSums(wmatrix)

  # Optimality function Value of the ideal alternative (S0)
  nideal <- numeric(m)
  for (j in 1:m) {
    nideal[j] <- max(nmatrix[, j])
  }
  
  S0 <- sum(nideal * weights)

  # Utility degree of alternatives (Ki)
  Ki <- Si / S0

  ranking <- rank(-Ki, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list( Si = Si, Ki = Ki, rank = ranking)

  return(results)
}
