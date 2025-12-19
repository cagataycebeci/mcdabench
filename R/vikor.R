vikor <- function(dmatrix, bcvec, weights, normethod = "maxmin", v = 0.5, tiesmethod="average") {
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
  
  # Normalize the decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = normethod)
  
  # Determine the best and worst values for each criterion
  f_star <- numeric(m)
  f_minus <- numeric(m)
  
  for (j in 1:m) {
    if (bcvec[j] == 1) {
      f_star[j] <- max(dmatrix[, j])
      f_minus[j] <- min(dmatrix[, j])
    } else if (bcvec[j] == -1) {
      f_star[j] <- min(dmatrix[, j])
      f_minus[j] <- max(dmatrix[, j])
    } else {
      stop("Invalid criteria type. Must be 1 or -1.")
    }
  }
  
  S <- numeric(n)
  R <- numeric(n)
  
  for (i in 1:n) {
    S[i] <- sum(weights * (f_star - dmatrix[i, ]) / (f_star - f_minus))
    R[i] <- max(weights * (f_star - dmatrix[i, ]) / (f_star - f_minus))
  }
  
  S_star <- min(S)
  S_minus <- max(S)
  R_star <- min(R)
  R_minus <- max(R)
  
  Q <- numeric(n)
  for (i in 1:n) {
    Q[i] <- v * (S[i] - S_star) / (S_minus - S_star) + 
      (1 - v) * (R[i] - R_star) / (R_minus - R_star)
  }
  
  aranks <- rank(Q, ties.method = tiesmethod)
  names(aranks) <- rownames(dmatrix)
  
  results <- list(
    fminus = f_minus,
    fstar = f_star,
    S = S,
    R = R,
    Q = Q,
    rank = aranks
  )
  
  return(results)
}
