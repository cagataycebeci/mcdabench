mairca <- function(dmatrix, bcvec, weights, tiesmethod="average") {
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
  
  # Normalized decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = "maxmin")

  # Calculate PAj values
  PA <- rep(1 / m, m)

  # Calculate TPij values
  tpmatrix <- matrix(PA * weights, nrow = n, ncol = m, byrow = TRUE)

  # Calculate TRij values
  trmatrix <- tpmatrix * nmatrix 

  # Calculate Gij values
  gmatrix <- pmax(tpmatrix - trmatrix, 0) 

  # Calclate Qi values
  Q <- rowSums(gmatrix)

  ranking <- rank(-Q, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)

  results <- list(Q = Q, rank=ranking)

  return(results)
}
