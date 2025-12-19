ocra <- function(dmatrix, bcvec, weights, tiesmethod = "average") {
  # Input Checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(1, -1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain 1 (benefit) or -1 (cost) values.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Determine Reference Values (Max for Benefit, Min for Cost)
  rvals <- numeric(m)
  for (j in seq_len(m)) {
    rvals[j] <- ifelse(bcvec[j] == 1, max(dmatrix[, j]), min(dmatrix[, j]))
  }

  # Prevent division errors
  rvals[rvals == 0] <- 1e-6

  # Compute Impact Values
  I <- rep(0, n)  # Cost impact
  O <- rep(0, n)  # Benefit impact
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      if (bcvec[j] == 1) {
        O[i] <- O[i] + weights[j] * ((rvals[j] - dmatrix[i, j]) / rvals[j])
      } else {
        I[i] <- I[i] + weights[j] * ((dmatrix[i, j] - rvals[j]) / rvals[j])
      }
    }
  }

  # Normalize I and O
  O <- O - min(O)
  I <- I - min(I)

  # Calculate Total Performance Score (P)
  P <- I + O - min(I + O)

  ranking <- rank(P, ties.method = tiesmethod)
  
  results <- list(
    reference_values = rvals,
    cost_impact_I = I,
    benefit_impact_O = O,
    total_performance_P = P,
    rank = ranking
  )

  return(results)
}

