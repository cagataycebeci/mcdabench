electre1 <- function(dmatrix, bcvec, weights, ct = 0.65, dt = 0.35, tiesmethod = "average") {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("dmatrix must be a matrix or data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must match number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec must contain only -1 or 1.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = "maxmin")

  concordance <- matrix(0, n, n)
  discordance <- matrix(0, n, n)
  outrank <- matrix(FALSE, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      c_sum <- 0
      d_max <- 0
      for (k in 1:m) {
        d <- nmatrix[i, k] - nmatrix[j, k]
        if (bcvec[k] == -1) d <- -d
        if (d >= 0) c_sum <- c_sum + weights[k]
        if (d < 0 && abs(d) > d_max) d_max <- abs(d)
      }
      concordance[i, j] <- c_sum / sum(weights)
      discordance[i, j] <- d_max / max(abs(nmatrix[i, ] - nmatrix[j, ]))
      outrank[i, j] <- (concordance[i, j] >= ct) && (discordance[i, j] <= dt)
    }
  }

  # Kernel (alternatives not strongly dominated)
  dominated <- apply(outrank, 2, any)
  kernel <- which(!dominated)
  names(kernel) <- rownames(dmatrix)[kernel]

  # Score = number of alternatives outranked
  outrank_score <- rowSums(outrank)
  ranking <- rank(-outrank_score, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)

  return(list(
    Concordance = concordance,
    Discordance = discordance,
    Outrank = outrank,
    Kernel = kernel,
    Outrank_Score = outrank_score,
    rank = ranking
  ))
}
