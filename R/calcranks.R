# Ranking calculaion methods
calcranks <- function(rankmatrix, rankmethod = "score") {
  n <- ncol(rankmatrix)  
  m <- nrow(rankmatrix)  

  # Initialize the A-score matrix (S_A in the paper)
  ascores <- matrix(0, nrow = n, ncol = n)
  colnames(ascores) <- colnames(rankmatrix)
  rownames(ascores) <- colnames(rankmatrix)

  # Populate the A-score matrix
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        ascores[i, j] <- 0 
        next
      }

      winsij <- 0 
      for (idx in 1:m) {
        if (rankmatrix[idx, i] < rankmatrix[idx, j]) { 
          winsij <- winsij + 1
        }
      }
      ascores[i, j] <- winsij
    }
  }

  # Calculate ranks based on the specified method
  # Rating skorlarını tutacak vektör
  rscores <- numeric(n)
  names(rscores) <- rownames(ascores) # Alternatif isimlerini ata

  if (rankmethod == "score") {
    for (i in 1:n) {
      rscores[i] <- sum(ascores[i, ])
      if (n > 0) {
        rscores[i] <- rscores[i] / n 
      } else {
        rscores[i] <- 0
      }
    }
  } else if (rankmethod == "neustadt") {
    for (i in 1:n) {
      sum_S_ij <- sum(ascores[i, ]) 
      sum_S_ji <- sum(ascores[, i])
      rscores[i] <- sum_S_ij - sum_S_ji
    }
  } else if (rankmethod == "buchholz") {
   for (i in 1:n) {
      sum_diff_sq <- 0
      for (j in 1:n) {
        if (i != j) { 
          sum_diff_sq <- sum_diff_sq + (ascores[i, j] - ascores[j, i])^2
        }
      }
      rscores[i] <- sum_diff_sq
    }
  } else {
    stop(paste0("Error: '", rankmethod, "' is an undefined or unsupported ranking method."))
  }
  
  rankscores <- if (rankmethod == "buchholz") rscores else -rscores

  score_ranks_min_tie <- rank(rankscores, ties.method = "min")
  names(score_ranks_min_tie) <- names(rankscores)

  average_ties_ranks <- rank(rankscores, ties.method = "average")
  names(average_ties_ranks) <- names(rankscores) 

  ordered_alternatives <- names(rscores)[order(rscores, decreasing = TRUE)]

  results <- list(
    rankmethod = rankmethod,
    rating_scores = rscores, 
    ordered_alternatives = ordered_alternatives,
    score_ranks = score_ranks_min_tie, 
    average_score_ranks = average_ties_ranks 
  )

  return(results)
}