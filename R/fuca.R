fuca <- function(dmatrix, bcvec, weights, tiesmethod="average") {
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
  if (sum(weights) != 1) {
    stop("The sum of the weights must be 1.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  rmatrix <- matrix(0, nrow = n, ncol = m)
  colnames(rmatrix) <- colnames(dmatrix)
  rownames(rmatrix) <- rownames(dmatrix)

  for (j in 1:m) {
    values <- dmatrix[, j]
    ranked_indices <- order(values, decreasing = (bcvec[j] == 1)) # Decreasing for benefit, increasing for cost
    ranks <- rank(values[ranked_indices], ties.method = tiesmethod)
    original_order_ranks <- numeric(n)
    original_order_ranks[ranked_indices] <- ranks
    rmatrix[, j] <- original_order_ranks
  }

  weighted_rank_scores <- rowSums(rmatrix * weights)
  
  ranking <- rank(-weighted_rank_scores, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(
    Rank_Matrix = rmatrix,
    Weighted_Rank_Scores = weighted_rank_scores,
    rank = ranking
  )
  return(results)
}
