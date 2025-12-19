mabac <- function(dmatrix, bcvec, weights, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("The sum of the weights must be approximately 1.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  # Normalization
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = "maxmin")

  # Weighted transformation
  wmatrix <- sweep(nmatrix + 1, 2, weights, "*")

  # Geometric mean (adjusted for equal influence of alternatives)
  G <- apply(wmatrix, 2, function(x) prod(x) ^ (1 / n))

  # Border approximation area (adjusted for numerical stability)
  Q <- sweep(wmatrix, 2, G, "-") / (1 + abs(G) + 1e-6)
  
  # Scores and ranking
  score <- rowSums(Q)
  ranking <- rank(-score, ties.method=tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  result <- list(G = G, Q = Q, score = score, rank=ranking)
  return(result)
}
