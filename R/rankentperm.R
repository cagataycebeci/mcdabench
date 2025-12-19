rankentperm <- function(rankmat, nperms = 1000, entropyopt = "shannon", padjmethod="none") {
  # Shannon Entropy
  shannon_entropy <- function(p) {
    return(-sum(p * log(p + 1e-10)))
  }

  # Jensen-Shannon Divergence
  jsd <- function(p1, p2) {
    if (length(p1) != length(p2)) stop("The compared vectors must have the same length.")

    x <- (p1 + p2) / 2
    jsdiv <- (sum(p1 * log((p1 + 1e-10) / (x + 1e-10))) +
              sum(p2 * log((p2 + 1e-10) / (x + 1e-10)))) / 2
    return(jsdiv)
  }

  n <- nrow(rankmat)

  entropy_matrix <- matrix(1, nrow = n, ncol = n, dimnames = list(rownames(rankmat), rownames(rankmat)))
  pval_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(rownames(rankmat), rownames(rankmat)))

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      rank1 <- as.numeric(rankmat[i, ])
      rank2 <- as.numeric(rankmat[j, ])

      normalized_rank1 <- (rank1 - min(rank1)) / (max(rank1) - min(rank1) + 1e-10)
      normalized_rank2 <- (rank2 - min(rank2)) / (max(rank2) - min(rank2) + 1e-10)

      # Entropy difference
      entropy_diff <- if (entropyopt == "jsd") {
        jsd(normalized_rank1, normalized_rank2)
      } else {
        abs(shannon_entropy(normalized_rank1) - shannon_entropy(normalized_rank2))
      }

      # Permutation test
      perm_diffs <- numeric(nperms)

      combined_ranks <- c(rank1, rank2)
      n_cols <- ncol(rankmat)

      for (k in 1:nperms) {
        perm_combined <- sample(combined_ranks)
        perm_rank1 <- perm_combined[1:n_cols]
        perm_rank2 <- perm_combined[(n_cols + 1):(2 * n_cols)]

        perm_normalized_rank1 <- (perm_rank1 - min(perm_rank1)) / (max(perm_rank1) - min(perm_rank1) + 1e-10)
        perm_normalized_rank2 <- (perm_rank2 - min(perm_rank2)) / (max(perm_rank2) - min(perm_rank2) + 1e-10)

        perm_diffs[k] <- if (entropyopt == "jsd") {
          jsd(perm_normalized_rank1, perm_normalized_rank2)
        } else {
          abs(shannon_entropy(perm_normalized_rank1) - shannon_entropy(perm_normalized_rank2))
        }
      }

      pval <- mean(perm_diffs >= entropy_diff)
      
      pval <- p.adjust(pval, method=padjmethod)

      entropy_matrix[i, j] <- entropy_diff
      entropy_matrix[j, i] <- entropy_diff
      pval_matrix[i, j] <- pval
      pval_matrix[j, i] <- pval
    }
  }

  displaymat <- matrix("", nrow = n, ncol = n, dimnames = list(rownames(rankmat), rownames(rankmat)))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        displaymat[i, j] <- format(round(pval_matrix[i, j], 2), nsmall = 2)
      } else if (i > j) {
        displaymat[i, j] <- format(round(entropy_matrix[i, j], 2), nsmall = 2)
      } else {
        displaymat[i, j] <- ""
      }
    }
  }

  displaymat <- as.data.frame(displaymat)

  return(list(
    nperms = nperms,
    entropyopt = entropyopt,
    entropy_matrix = entropy_matrix,
    p_adjustment = padjmethod,
    pval_matrix = pval_matrix,
    displaymat = displaymat
  ))
}
