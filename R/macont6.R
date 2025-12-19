# MACONT 6
macont6 <- function(dmatrix, bcvec, weights, p=0.5, q=0.5, delta=0.5, theta=0.5, tiesmethod="average") {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }

  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
      stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }

  m <- nrow(dmatrix)
  n <- ncol(dmatrix)

  nmatrix1 <- calcnormal(dmatrix, bcvec, type="ratio")
  nmatrix2 <- calcnormal(dmatrix, bcvec, type="maxmin")
  nmatrix3 <- calcnormal(dmatrix, bcvec, type="vector")

  # Calculate balanced calcnormald values
  xbmatrix <- p * nmatrix1 + q * nmatrix2 + (1 - p - q) * nmatrix3

  # Compute S1(a_i) and S2(a_i) values.
  S1 <- numeric(m)
  S2 <- numeric(m)

  # Calculate rho_i and zeta_i values
  rho <- numeric(m)
  zeta <- numeric(m)

  for (i in 1:m) {
    rho_sum_sq <- 0
    zeta_prod_num <- 1
    zeta_prod_den <- 1

    for (j in 1:n) {
      rho_sum_sq <- rho_sum_sq + (weights[j] * (xbmatrix[i, j] - mean(xbmatrix[, j])))^2

      if (xbmatrix[i, j] >= mean(xbmatrix[, j])) {
        zeta_prod_num <- zeta_prod_num * (xbmatrix[i, j] - mean(xbmatrix[, j]))^weights[j]
      }
      if (xbmatrix[i, j] <= mean(xbmatrix[, j])) {
        zeta_prod_den <- zeta_prod_den * (mean(xbmatrix[, j]) - xbmatrix[i, j])^weights[j]
      }
    }
    rho[i] <- sqrt(rho_sum_sq)
    if (zeta_prod_den != 0) {
      zeta[i] <- zeta_prod_num / zeta_prod_den
    } else {
      zeta[i] <- 0 # Avoid division by zero
    }

    S1[i] <- delta * rho[i] / sqrt(sum(rho^2)) + (1 - delta) * zeta[i] / sqrt(sum(zeta^2))
    S2[i] <- theta * max(weights * (xbmatrix[i, ] - colMeans(xbmatrix))) +
             (1 - theta) * min(weights * (xbmatrix[i, ] - colMeans(xbmatrix)))
  }

  # Compute final comprehensive score
  S <- (S1 + S2) / (2 + sqrt(sum(S1^2)))

  # Rank alternatives based on comprehensive scores
  # ranking <- order(S, decreasing = TRUE)
  ranking <- rank(-S, ties.method=tiesmethod)
  
  results <- list(
    balance_matrix = xbmatrix,
    S1_values = S1,
    S2_values = S2,
    comprehensive_scores = S,
    rank = ranking
  )
  return(results)
}
