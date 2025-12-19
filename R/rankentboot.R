rankentboot <- function(rankmat, nboot = 1000, 
  entropyopt = "shannon", padjmethod = "none") {

  if (!(padjmethod %in% p.adjust.methods)) {
     stop("Invalid padjmethod. Choose from: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'.")
  }

  n <- nrow(rankmat)
  result_matrix <- matrix("", nrow = n, ncol = n) 
  rownames(result_matrix) <- rownames(rankmat)
  colnames(result_matrix) <- rownames(rankmat)

  # Shannon entropy calculation
  shannon_entropy <- function(ranks) {
    normalized_ranks <- (ranks - min(ranks)) / (max(ranks) - min(ranks) + 1e-10)
    return(-sum(normalized_ranks * log(normalized_ranks + 1e-10)))
  }

  # Jensen-Shannon divergence calculation
  jsd_divergence <- function(x, y) {
    px <- (x - min(x)) / (max(x) - min(x) + 1e-10)
    py <- (y - min(y)) / (max(y) - min(y) + 1e-10)
    m <- (px + py) / 2

    kl_x <- sum(px * log(px / (m + 1e-10) + 1e-10))
    kl_y <- sum(py * log(py / (m + 1e-10) + 1e-10))

    return(0.5 * (kl_x + kl_y))
  }

  # Tsallis entropy function
  tsallis_entropy <- function(ranks, q) {
    if (q == 1) {
      return(-sum(ranks * log(ranks)))
    } else {
      return((1 - sum(ranks^q)) / (q - 1))
    }
  }

  raw_p_values <- numeric()
  indices <- list()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      x <- rankmat[i, ]
      y <- rankmat[j, ]

      if (entropyopt == "shannon") {
        entropy_x_orig <- shannon_entropy(x)
        entropy_y_orig <- shannon_entropy(y)
        obs_stat <- abs(entropy_x_orig - entropy_y_orig)
      } else if (entropyopt == "jsd") {
        obs_stat <- jsd_divergence(x, y)
      } else if (entropyopt == "tsallis") {
        entropy_x_orig <- tsallis_entropy(x, q=2)
        entropy_y_orig <- tsallis_entropy(y, q=2)
        obs_stat <- abs(entropy_x_orig - entropy_y_orig)
      } else {
        stop("Invalid entropyopt. Use 'shannon' 'tsallis' or 'jsd'.")
      }

      # Compute similarity score (inverse entropy difference)
      similarity_score <- format(round(1 - obs_stat, 3), nsmall = 3)
      result_matrix[j, i] <- similarity_score 

      # Bootstrap process for p-value calculation
      boot_stats <- numeric(nboot)
      n_combined <- length(c(x, y))

      for (b in 1:nboot) {
        combined_boot_indices <- sample(1:n_combined, size = n_combined, replace = TRUE)
        boot_x <- c(x, y)[combined_boot_indices[1:length(x)]]
        boot_y <- c(x, y)[combined_boot_indices[(length(x) + 1):n_combined]]

        if (entropyopt == "shannon") {
          boot_entropy_x <- shannon_entropy(boot_x)
          boot_entropy_y <- shannon_entropy(boot_y)
          boot_stats[b] <- abs(boot_entropy_x - boot_entropy_y)
        } else {
          boot_stats[b] <- jsd_divergence(boot_x, boot_y)
        }
      }

      # Store raw p-values for adjustment
      raw_p_values <- c(raw_p_values, mean(boot_stats >= obs_stat))
      indices <- c(indices, list(c(i, j)))
    }
  }

   adj_p_values <- p.adjust(raw_p_values, method = padjmethod)

   for (idx in seq_along(adj_p_values)) {
    result_matrix[indices[[idx]][1], indices[[idx]][2]] <- format(round(adj_p_values[idx], 3), nsmall = 3)
  }

  display_matrix <- as.data.frame(result_matrix)
  results <- list(entropy=entropyopt, p_adjustment=padjmethod, nboot=nboot, displaymat = display_matrix)
  return(results)
}
