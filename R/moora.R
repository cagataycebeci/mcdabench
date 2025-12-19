moora <- function(dmatrix, bcvec, weights, tiesmethod="average") {
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

  # Normalize Decision Matrix
  vmatrix <- matrix(0, nrow = n, ncol = m)
  for (j in 1:m) {
    denominator <- sqrt(sum(dmatrix[, j]^2))
    if (denominator == 0) {
      vmatrix[, j] <- 0
    } else {
      vmatrix[, j] <- dmatrix[, j] / denominator
    }
  }

  # Weighted Normalized Decision Matrix
  wmatrix <- vmatrix * weights

  # Optimization Values
  optvals <- numeric(n)
  for (i in 1:n) {
    bsum <- sum(wmatrix[i, which(bcvec == 1)])  # Benefit sum
    csum <- sum(wmatrix[i, which(bcvec == -1)]) # Cost sum
    optvals[i] <- bsum - csum
  }

  ranking <- rank(-optvals, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(
    Normalized_Matrix = vmatrix,
    Weighted_Normalized_Matrix = wmatrix,
    Optimization_Values = optvals,
    rank = ranking
  )
  return(results)
}