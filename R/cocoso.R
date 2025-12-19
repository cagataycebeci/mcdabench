cocoso <- function(dmatrix, bcvec, weights, lambda = 0.5, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec (benefit/cost vector) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (lambda < 0 || lambda > 1) {
    stop("The lambda must be between 0 and 1.")
  }
  
  # Normalize the decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type="linear")
  
  wmatrix <- sweep(nmatrix, 2, weights, "*")
  
  # Calculate Si and Pi values
  Si <- rowSums(wmatrix)
  Pi <- apply(wmatrix, 1, function(row) prod(row^weights))
 
  # Calculate the appraisal score strategies
  # kia <- (Pi + Si) / (sum(Pi + Si))
  kia <- Si
  kib <- Si / min(Si) + Pi / min(Pi)
  kic <- (lambda * Si + (1 - lambda) * Pi) / (lambda * max(Si) + (1 - lambda) * max(Pi))
  
  # Calculate the performance scores of the alternatives
  Ki <- (kia * kib * kic)^(1 / 3) + (1 / 3) * (kia + kib + kic)
  
  ranking <- rank(-Ki, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(S_i = Si, P_i = Pi, k_ia = kia, k_ib = kib, k_ic = kic, k_i = Ki, rank = ranking)
  
  return(results)
}

