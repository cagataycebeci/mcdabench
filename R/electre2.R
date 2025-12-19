electre2 <- function(dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, 
            st = 0.65, wt = 0.5, tiesmethod = "average") {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("dmatrix must be a matrix or data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must match the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec must contain only -1 or 1 values.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  nmatrix <- calcnormal(dmatrix, bcvec = bcvec, type = "maxmin")

  if (is.null(p)) p <- rep(list(list(p = max(abs(apply(nmatrix, 2, diff))))), m)
  if (is.null(q)) q <- rep(list(list(q = p[[1]]$p / 2)), m)
  if (is.null(v)) v <- rep(list(list(v = p[[1]]$p * 1.5)), m)

  concordance_matrix <- matrix(0, n, n)
  discordance_matrix <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      concordance_sum <- 0
      discordance_max <- 0
      for (k in 1:m) {
        d <- nmatrix[i, k] - nmatrix[j, k]
        if (bcvec[k] == -1) d <- -d

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
      discordance_matrix[i, j] <- ifelse(discordance_max > 0, discordance_max / max(abs(nmatrix[i, ] - nmatrix[j, ])), 0)
    }
  }

  credibility_matrix <- concordance_matrix * (1 - discordance_matrix)
  strong_rel <- (credibility_matrix >=st)
  weak_rel   <- (credibility_matrix >= wt) & !strong_rel

  strong_out <- rowSums(strong_rel)
  strong_in  <- colSums(strong_rel)
  weak_out   <- rowSums(weak_rel)
  weak_in    <- colSums(weak_rel)

  upward_score <- strong_out + 0.5 * weak_out
  downward_score <- strong_in + 0.5 * weak_in

  # Apply tiesmethod in final ranking step
  ranking_up <- rank(-upward_score, ties.method = tiesmethod)
  ranking_down <- rank(downward_score, ties.method = tiesmethod)
  average_rank <- (ranking_up + ranking_down) / 2
  ranking <- rank(average_rank, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)

  return(list(
    Concordance = concordance_matrix,
    Discordance = discordance_matrix,
    Credibility = credibility_matrix,
    Strong = strong_rel,
    Weak = weak_rel,
    Rank_Upward = ranking_up,
    Rank_Downward = ranking_down,
    rank = ranking
  ))
}
