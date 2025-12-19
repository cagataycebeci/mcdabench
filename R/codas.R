codas <- function(dmatrix, bcvec, weights, thr = 0.1, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matirx (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (thr < 0 || thr > 1) {
    stop("The threshold (thr) must be between 0 and 1.")
  }
  
  # Normaliz the decision matrix
  nmatrix <- calcnormal(X=dmatrix, bcvec=bcvec, type="ratio")
  
  wmatrix <- sweep(nmatrix, 2, weights, "*")

  # Ideal solution
  niv <- numeric(ncol(wmatrix))
  for (j in 1:ncol(wmatrix)) {
    if (bcvec[j] == 1) { 
      niv[j] <- min(wmatrix[, j])
    } else { 
      niv[j] <- max(wmatrix[, j])
    }
  }
  
  # Calculate Euclidean and Taxicab distances
  edist <- apply(wmatrix, 1, function(x) sqrt(sum((x - niv)^2)))
  tdist <- apply(wmatrix, 1, function(x) sum(abs(x - niv)))

  # Calculate appraisal score (Pi)
  pii <- edist + (1 - thr) * (tdist)

  ranking <- rank(-pii, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(niv = niv, euc_dist = edist, tax_dist = tdist, pi_i = pii, rank = ranking)
  return(results)
}
