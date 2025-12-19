# Rank Aggregation
rankaggregate <- function(rankmat, topk=3, damp=0.5, niters=100, tol=1e-4, tiesmethod="average") {

  # Top-K Counts Aggregation 
  topkcounts <- function(rankmat, topk=1, ties="average") {
      topkcounts <- sapply(colnames(rankmat), function(x) sum(rankmat[, x] <= topk + 0.5))
      names(topkcounts) <- colnames(rankmat)
      topk_ranking <- rank(-topkcounts, ties.method=ties)
      results <- list(topk=topk, topk_counts=topkcounts, topk_ranking=topk_ranking)
      return(results)
  }

  # Rank Sum Aggregation
  ranksum <- function(rankmat, ties="average"){
     ranksums <- colSums(rankmat)
     names(ranksums) <- colnames(rankmat)
     ranksum_ranking <- rank(ranksums, ties.method = ties)
     results <- list(rank_sums=ranksums, ranksum_ranking=ranksum_ranking)
     return(results)
  }
  
  # Copeland Rank Aggregation
  copeland <- function(rankmat, ties = "average") {  
    n <- nrow(rankmat)
    m <- ncol(rankmat) 
    cnames <- colnames(rankmat)
    
    copeland_scores <- numeric(m)
    names(copeland_scores) <- cnames
    
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
        wins_i <- sum(rankmat[, i] < rankmat[, j])
        ties_ij <- sum(rankmat[, i] == rankmat[, j])
        ratio_i <- (wins_i + 0.5 * ties_ij) / n
        
        wins_j <- sum(rankmat[, j] < rankmat[, i])
        ratio_j <- (wins_j + 0.5 * ties_ij) / n
        
        if (ratio_i > ratio_j) {
          copeland_scores[i] <- copeland_scores[i] + 1
        } else if (ratio_i < ratio_j) {
           copeland_scores[j] <- copeland_scores[j] + 1
        } else {
          copeland_scores[i] <- copeland_scores[i] + 0.5
          copeland_scores[j] <- copeland_scores[j] + 0.5
        }
      }
    }
    
    copeland_ranking <- rank(-copeland_scores, ties.method = ties)
    results <- list(copeland_scores = copeland_scores, copeland_ranking = copeland_ranking)
    return(results)
  }

  # Median Rank Aggregation
  rankmedian <- function(rankmat, ties="average") {
    medians <- apply(rankmat, 2, median)
    median_ranking <- rank(medians, ties.method = ties)
    results <- list(medians = medians, median_ranking = median_ranking)
    return(results)
  }

  # Recursive function to generate all permutations of a vector
  permute <- function(vec) {
    if (length(vec) == 1) {
      return(list(vec))
    } else {
      result <- list()
      # Loop over each element in vec
      for (i in seq_along(vec)) {
        # Remove i-th element to form the remaining vector
        remain <- vec[-i]
        # Recursively generate permutations of the remaining elements
        subperms <- permute(remain)
        # Append the removed element to the beginning of each sub-permutation
        for (sub in subperms) {
          result[[length(result) + 1]] <- c(vec[i], sub)
        }
      }
      return(result)
    }
  }
   
  kemenyyoung2 <- function(rankmat, ties="average") {
     m <- ncol(rankmat)
    candidates <- colnames(rankmat)
    if (is.null(candidates)) {
      candidates <- paste0("Alt", seq_len(m))
    }
    
    # Build the pairwise comparison matrix 'd'.
    # d[i, j] counts the number of voters that prefer candidate i over candidate j.
    d <- matrix(0, nrow = m, ncol = m, dimnames = list(candidates, candidates))
    for (i in 1:m) {
      for (j in 1:m) {
        if (i != j) {
          d[i, j] <- sum(rankmat[, candidates[i]] < rankmat[, candidates[j]])
        }
      }
    }
    
    # Generate all possible permutations (rankings) of the candidates.
    all_perms <- permute(candidates)
    
    best_score <- -Inf
    best_perms <- list()  # Store all permutations that achieve the best (highest) score.
    
    # For each permutation, compute its Kemeny score.
    # For permutation 'pi', the Kemeny score is the sum of d[pi[i], pi[j]]
    # for every pair (i, j) with i < j.
    for (perm in all_perms) {
      score <- 0
      for (i in 1:(m - 1)) {
        for (j in (i + 1):m) {
          score <- score + d[perm[i], perm[j]]
        }
      }
      # Update the best score and track all permutations that match it.
      if (score > best_score) {
        best_score <- score
        best_perms <- list(perm)
      } else if (score == best_score) {
        best_perms[[length(best_perms) + 1]] <- perm
      }
    }
    
    # Sum the positions for each candidate over all best permutations.
    # 'ranking_sum' collects each candidate's position (e.g., first position, second, etc.)
    # across the best permutations.
    ranking_sum <- setNames(rep(0, m), candidates)
    for (perm in best_perms) {
      for (pos in seq_along(perm)) {
        ranking_sum[perm[pos]] <- ranking_sum[perm[pos]] + pos
      }
    }
    
    # Compute the average position for each candidate by dividing by the number of best permutations.
    avg_positions <- ranking_sum / length(best_perms)
    avg_ranking <- rank(avg_positions, ties.method = ties)
    
    results <- list(
      kemeny_ranking = avg_ranking,
      average_positions = avg_positions,
      kemeny_score = best_score,
      pairwise_matrix = d,
      best_perms = best_perms
    )
    
    return(results)
  }
  
  # Kemeny & Young Rank Aggregation II
  # (based on Kendall Tau distance)
  kemenyyoung <- function(rankmat, ties="average") {
      m <- ncol(rankmat)
      cnames <- colnames(rankmat)
      n <- nrow(rankmat)
  
      preference_matrix <- matrix(0, nrow = m, ncol = m, dimnames = list(cnames, cnames))
  
      for (i in 1:m) {
        for (j in 1:m) {
          if (i != j) {
            wins <- sum(rankmat[, i] < rankmat[, j])
            losses <- sum(rankmat[, i] > rankmat[, j])
            preference_matrix[i, j] <- wins - losses
          }
        }
      }
  
      kemeny_scores <- numeric(m)
      names(kemeny_scores) <- cnames
      for (i in 1:m) {
        for (j in 1:m) {
          if (i != j) {
            kemeny_scores[i] <- kemeny_scores[i] + sign(preference_matrix[i, j])
          }
        }
      }
  
      kemeny_ranking <- rank(-kemeny_scores, ties.method = ties)
  
      results <- list(kemeny_scores = kemeny_scores, kemeny_ranking = kemeny_ranking)
      return(results)
  }

  # Markov Chain Rank Aggreation
  markov <- function(rankmat, damping = 0.85, tol = 1e-9, max_iter = 1000, ties = "average") {
    m <- ncol(rankmat)
    n <- nrow(rankmat)
    cnames <- colnames(rankmat)
    wins <- matrix(0, nrow = m, ncol = m, dimnames = list(cnames, cnames))
  
    for (k in 1:n) {
      ranks <- rank(rankmat[k, ], ties.method = ties)
      for (i in 1:m) {
        for (j in 1:m) {
          if (i != j) {
            if (ranks[i] < ranks[j]) {
              wins[i, j] <- wins[i, j] + 1
            }
            # Optional: Handling ties (example - you might need a different approach)
            # else if (ranks[i] == ranks[j]) {
            #   # Could add a smaller value or do nothing
            # }
          }
        }
      }
    }
  
    T <- matrix(0, nrow = m, ncol = m, dimnames = list(cnames, cnames))
    for (j in 1:m) {
      col_sum <- sum(wins[, j])
      if (col_sum > 0) {
        T[, j] <- wins[, j] / col_sum
      } else {
        T[, j] <- rep(1/m, m)
      }
    }
  
    v <- rep(1/m, m) # Initial uniform vector
    names(v) <- cnames
  
    # Optional: Initialize with prior knowledge (if available)
    # initial_v <- ...
    # if (!is.null(initial_v) && length(initial_v) == m && all(names(initial_v) == cnames)) {
    #   v <- initial_v
    # }
  
    for (iter in 1:max_iter) {
      new_v <- damping * (T %*% v) + (1 - damping) / m
      # Optional: Different convergence criterion (relative change)
      # if (sum(abs((new_v - v) / v)) < tol) break
      if (sum(abs(new_v - v)) < tol) break
      v <- new_v
    }
  
    markov_scores <- as.vector(v)
    names(markov_scores) <- cnames
  
    # Optional: Normalize scores to a specific range (e.g., 0 to 1)
    # markov_scores_normalized <- (markov_scores - min(markov_scores)) / (max(markov_scores) - min(markov_scores))
    # markov_ranking <- rank(-markov_scores_normalized, ties.method = ties)
  
    markov_ranking <- rank(-markov_scores, ties.method = ties)
  
    return(list(iterations = iter,
                wins_matrix = wins,
                transition_matrix = T,
                markov_scores = markov_scores,
                markov_ranking = markov_ranking))
    }
    
    # Borda Count ranking
    bordacount <- function(rankmat, ties="average") {
       n <- nrow(rankmat)
       m <- ncol(rankmat)
   
       borda_scores <- numeric(m) 
       names(borda_scores) <- colnames(rankmat)  
  
       for (i in 1:n) {  
          borda_scores <- borda_scores + (n - rankmat[i, ])
       }
  
       borda_ranking <- rank(-borda_scores, ties.method = ties)  
  
       return(list(borda_scores = borda_scores, borda_ranking = borda_ranking))
    }

    # Outranking table
    outranking_table <- function(prefmatrix) {
      outranking <- function(ranks) {
        if (is.null(names(ranks))) {
           names(ranks) <- paste0("ALT_", seq_along(ranks))
        }
        sranks <- sort(ranks)
        hierarchy <- ""
        for (i in seq_along(sranks)) {
           alt <- names(sranks)[i]
           hierarchy <- paste0(hierarchy, alt)
           if (i < length(sranks)) {
              if (sranks[i] == sranks[i + 1]) {
                 hierarchy <- paste0(hierarchy, " = ")
              } else {
                hierarchy <- paste0(hierarchy, " > ")
              }
           }
        }
        return(hierarchy)
     }
  
     outranks <- apply(prefmatrix, 1, outranking)
     otable <- data.frame(
        Method = rownames(prefmatrix),
        Outranking = outranks,
        row.names = NULL,
        stringsAsFactors = FALSE
     )
     return(otable)
   }
  
   restopk <- topkcounts(rankmat, topk=topk, ties=tiesmethod)
   resranksum <- ranksum(rankmat, ties=tiesmethod)
   resmedian <- rankmedian(rankmat, ties=tiesmethod)
   rescopeland <- copeland(rankmat, ties=tiesmethod)
 # reskemeny2 <- kemenyyoung2(rankmat, ties=tiesmethod)
   reskemeny <- kemenyyoung(rankmat, ties=tiesmethod)
   resmarkov <- markov(rankmat, damping = damp, tol=tol, max_iter=niters, ties=tiesmethod)
   resborda <- bordacount(rankmat, ties=tiesmethod)
 # ressrd <- ranksrd(rankmat, refopt = "mean")
   prefmat <- rbind(restopk$topk_ranking, resranksum$ranksum_ranking, resmedian$median_ranking, resborda$borda_ranking, 
      rescopeland$copeland_ranking, reskemeny$kemeny_ranking, resmarkov$markov_ranking) 
   rownames(prefmat) <- c(paste0("TOPK",topk), "RANKSUM", "MEDIAN", "BORDACNT", "COPELAND", "KEMYNG", "MARKOV") 
   preftable <- outranking_table(prefmat)
  
   results <- list(TOPK=restopk, RANKSUM=resranksum, MEDIAN=resmedian, BORDA_COUNT=resborda, 
      COPELAND=rescopeland, KEMENY_YOUNG=reskemeny, MARKOV=resmarkov,  
      preference_ranking=prefmat, preference_table=preftable) 

   return(results)
}
