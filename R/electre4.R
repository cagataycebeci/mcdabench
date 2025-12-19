electre4 <- function(dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, thr = 0.5, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec (benefit/cost vector) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Normalized decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type="maxmin")

  # Default pseudo-criterion thresholds
  if (is.null(p)) {
    p <- rep(list(list(p = max(abs(apply(nmatrix, 2, diff))))), m)
  }
  if (is.null(q)) {
    q <- rep(list(list(q = p[[1]]$p / 2)), m)
  }
  if (is.null(v)) {
    v <- rep(list(list(v = p[[1]]$p * 1.5)), m)
  }

  # Concordance & Discordance Matrices
  concordance_matrix <- matrix(0, nrow = n, ncol = n)
  discordance_matrix <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        concordance_sum <- 0
        discordance_max <- 0
        for (k in 1:m) {
          d <- nmatrix[i, k] - nmatrix[j, k]
          if (bcvec[k] == -1) {
            d <- -d
          }

          # Pseudo-criterion evaluation
          if (d >= q[[k]]$q) {
            concordance_sum <- concordance_sum + weights[k]
          }
          if (d < -v[[k]]$v) { 
            discordance_val <- abs(d)
            if (discordance_val > discordance_max) {
              discordance_max <- discordance_val
            }
          }
        }
        concordance_matrix[i, j] <- concordance_sum / sum(weights) 
        #discordance_matrix[i, j] <- ifelse(discordance_max > 0, discordance_max / max(abs(nmatrix[, k])), 0)
 	discordance_matrix[i, j] <- ifelse(discordance_max > 0, discordance_max / max(abs(nmatrix[i, ] - nmatrix[j, ])), 0)

      } else {
        concordance_matrix[i, j] <- 1
        discordance_matrix[i, j] <- 0
      }
    }
  }

  # Credibility matrix
  credibility_matrix <- concordance_matrix * (1 - discordance_matrix)

  # Elimination process
  dominance_matrix <- credibility_matrix > thr
  scores <- rowMeans(credibility_matrix) 

  ranking <- rank(-scores, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)

  results <- list(
     Concordance = concordance_matrix,
     Discordance = discordance_matrix,
     Credibility = credibility_matrix,
     Net_Flow = scores,
     rank = ranking
  )

  return(results)
}
