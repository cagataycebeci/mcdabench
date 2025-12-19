electre3 <- function(dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, thr = 0.5, tiesmethod = "average") {
  # Input validation
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("dmatrix must be a matrix or data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec must contain -1 (cost) or 1 (benefit) values and match the number of criteria.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Normalize the decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = "maxmin")

  # Generate default thresholds if not provided
  if (is.null(p)) p <- rep(list(list(p = max(abs(apply(nmatrix, 2, diff))))), m)
  if (is.null(q)) q <- rep(list(list(q = p[[1]]$p / 2)), m)
  if (is.null(v)) v <- rep(list(list(v = p[[1]]$p * 1.5)), m)

  credibility_matrix <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        credibility_matrix[i, j] <- 1
        next
      }

      concordance_sum <- 0
      partial_concordance <- numeric(m)
      discordance_vals <- numeric(m)

      for (k in 1:m) {
        d <- nmatrix[i, k] - nmatrix[j, k]
        if (bcvec[k] == -1) d <- -d

        dq <- q[[k]]$q
        dp <- p[[k]]$p
        dv <- v[[k]]$v

        # Fuzzy concordance computation
        if (d <= dq) {
          ckj <- 0
        } else if (d <= dp) {
          ckj <- (d - dq) / (dp - dq)
        } else {
          ckj <- 1
        }
        partial_concordance[k] <- ckj
        concordance_sum <- concordance_sum + ckj * weights[k]

        # Discordance value based on veto threshold
        if (d >= -dv) {
          djk <- 0
        } else {
          djk <- abs(d) / (abs(d) + 1e-6)  # Prevent division by zero
        }
        discordance_vals[k] <- djk
      }

      # Aggregate concordance index
      Cij <- concordance_sum / sum(weights)

      # Apply credibility weakening
      weakening <- 1
      for (k in 1:m) {
        if (discordance_vals[k] > Cij) {
          weakening <- weakening * (1 - (discordance_vals[k] - Cij) / (1 - Cij))
        }
      }

      # Final credibility score
      sigma_ij <- Cij * weakening
      credibility_matrix[i, j] <- sigma_ij
    }
  }

  # Aggregate credibility values and derive rankings
  scores <- rowMeans(credibility_matrix)
  ranking <- rank(-scores, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)

  return(list(
    Credibility = credibility_matrix,
    Net_Flow = scores,
    rank = ranking
  ))
}
