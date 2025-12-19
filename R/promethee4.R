promethee4 <- function(dmatrix, bcvec, weights, normethod = NULL, prefuncs = NULL, thr = NULL, alpha = 0.2, tiesmethod = "average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec (benefit/cost vector) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("The sum of the weights must be approximately 1.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  if (is.null(normethod)) {
    nmatrix <- dmatrix
  } else {
    nmatrix <- calcnormal(dmatrix, bcvec, type = normethod)
  }

  # Default preference functions (linear)
  if (is.null(prefuncs)) {
    prefuncs <- rep("linear", m)
  } else if (is.character(prefuncs)) {
    if (length(prefuncs) == 1) {
      prefuncs <- rep(prefuncs, m)
    } else if (length(prefuncs) != m) {
      stop("The length of prefuncs must be 1 or equal to the number of criteria.")
    }
  } else {
    stop("Invalid prefuncs! Must be a character vector or a single character value.")
  }

  # Handle thresholds
  if (is.null(thr)) {
    # thr <- as.list(rep(list(p = sd(dmatrix)), m))
    thr <- as.list(rep(list(p = max(abs(apply(dmatrix, 2, diff)))), m)) # Default threshold per criterion
    warning("No threshold provided. Using default threshold based on the maximum absolute difference for each criterion.")
  } else if (is.list(thr)) {
    if (length(thr) == 1 && m > 1) {
      thr <- rep(thr, m) # Replicate single threshold list for all criteria
    } else if (length(thr) != m) {
      stop("The length of threshold list must be 1 or equal to the number of criteria.")
    }
  } else {
    stop("Invalid threshold! It must be a list.")
  }

  # Calculate preference indices
  preference_indices <- array(0, dim = c(n, n, m)) # [alt i, alt j, criterion k]
  for (k in 1:m) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          d <- dmatrix[i, k] - dmatrix[j, k]
          # If the criterion requires to be minimized, reverse the sign of the difference
          if (bcvec[k] == -1) {
            d <- -d
          }

          # Preference functions
          pref_type <- prefuncs[k]
          thresh <- thr[[k]]

          if (pref_type == "linear") {
            if (is.list(thresh) && !is.null(thresh$p)) {
              p <- thresh$p
            } else if (is.numeric(thresh) && length(thresh) == 1) {
              p <- thresh
            } else {
              stop(paste("Missing or invalid 'p' threshold for linear preference function for criterion", k))
            }
            preference_indices[i, j, k] <- ifelse(d <= 0, 0, ifelse(d > p, 1, d / p))
          } else if (pref_type == "usual") {
            preference_indices[i, j, k] <- ifelse(d <= 0, 0, 1)
          } else if (pref_type == "quasi-criterion") {
            if (is.list(thresh) && !is.null(thresh$s)) {
              s <- thresh$s
            } else {
              stop(paste("Missing 's' threshold for quasi-criterion preference function for criterion", k))
            }
            preference_indices[i, j, k] <- ifelse(d <= s, 0, 1)
          } else if (pref_type == "v-shape") {
            if (is.list(thresh) && !is.null(thresh$p)) {
              p <- thresh$p
            } else if (is.numeric(thresh) && length(thresh) == 1) {
              p <- thresh
            } else {
              stop(paste("Missing or invalid 'p' threshold for v-shape preference function for criterion", k))
            }
            preference_indices[i, j, k] <- ifelse(d <= 0, 0, d / p)
          } else if (pref_type == "level") {
            if (is.list(thresh) && !is.null(thresh$p) && !is.null(thresh$q)) {
              p <- thresh$p
              q <- thresh$q
            } else if (is.numeric(thresh) && length(thresh) == 2) {
              p <- thresh[1]
              q <- thresh[2]
            }
            else {
              stop(paste("Missing 'p' or 'q' thresholds for level preference function for criterion", k))
            }
            preference_indices[i, j, k] <- ifelse(d <= q, 0, ifelse(d > p, 1, 0.5))
          } else {
            stop(paste("Invalid preference function:", pref_type, "for criterion", k))
          }
        }
      }
    }
  }

  # Weighted preference indices
  weighted_preference_indices <- apply(preference_indices, c(1, 2), function(x) sum(x * weights))

  # Calculate positive and negative outranking flows
  phi_plus <- rowSums(weighted_preference_indices) / (n - 1)
  phi_minus <- colSums(weighted_preference_indices) / (n - 1)

  # Calculate net outranking flow (PROMETHEE II)
  phi_net <- phi_plus - phi_minus

  # Calculate the average outranking flow
  phi_avg <- mean(phi_net)

  # Calculate the PROMETHEE IV ranking
  complete_ranking <- rank(-(phi_net - alpha * phi_avg), ties.method = tiesmethod)
  names(complete_ranking) <- rownames(dmatrix)
  
  results <- list(
    Phi_plus = phi_plus,
    Phi_minus = phi_minus,
    Phi_net = phi_net,
    Phi_avg = phi_avg,
    rank = complete_ranking
  )
  return(results)
}
