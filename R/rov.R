# Rank of Value (ROV) method
#
rov <- function(dmatrix, bcvec, weights, normethod="maxmin", tiesmethod="average") {
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
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("The sum of the weights must be approximately 1.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Normalize the decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type=normethod)
  
  ui <- rowSums(nmatrix * weights)

  ranking <- rank(-ui, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(u_i = ui, rank = ranking)

  return(results)
}
