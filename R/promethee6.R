promethee6 <- function(dmatrix, bcvec, weights, normethod = NULL, prefuncs = NULL, thr = NULL, varmethod = "abs_sum", p = NULL, q = NULL, tiesmethod = "average") {
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
  
  # Default Preference Functions
  if (is.null(prefuncs)) {
    prefuncs <- rep("linear", m)
  } else if (is.character(prefuncs)) {
    if (length(prefuncs) == 1) {
      prefuncs <- rep(prefuncs, m)
    } else if (length(prefuncs) != m) {
      stop("The length of prefuncs vector must be 1 or equal to the number of criteria.")
    }
  } else {
    stop("Invalid prefuncs! Must be a character vector or a single character value.")
  }

  # Handle Thresholds
  if (is.null(thr)) {
    warning("No threshold provided. PROMETHEE method results may be affected.")
    # thr <- as.list(rep(list(p = sd(dmatrix)), m))
    thr <- as.list(rep(list(p = max(abs(apply(dmatrix, 2, diff)))), m)) # A less 'surprising' default
  } else if (is.list(thr)) {
    if (length(thr) == 1 && m > 1) {
      thr <- rep(thr, m)
    } else if (length(thr) != m) {
      stop("The length of thr list must be 1 or equal to the number of criteria.")
    }
  } else {
    stop("Invalid thr! Must be a list.")
  }

  # Calculate Preference Indices
  preference_indices <- array(0, dim = c(n, n, m))
  for (k in 1:m) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          d <- dmatrix[i, k] - dmatrix[j, k]
          if (bcvec[k] == -1) {
            d <- -d
          }

          pref_type <- prefuncs[k]
          thresh <- thr[[k]]

          preference_indices[i, j, k] <- switch(pref_type,
            "linear" = {
              p_val <- if (is.list(thresh) && !is.null(thresh$p)) thresh$p else if (is.numeric(thresh) && length(thresh) == 1) thresh else stop(paste("Missing or invalid 'p' threshold for linear preference function (criterion", k, ")"))
              ifelse(d <= 0, 0, ifelse(d > p_val, 1, d / p_val))
            },
            "usual" = ifelse(d <= 0, 0, 1),
            "quasi-criterion" = {
              s_val <- if (is.list(thresh) && !is.null(thresh$s)) thresh$s else stop(paste("Missing 's' threshold for quasi-criterion preference function (criterion", k, ")"))
              ifelse(d <= s_val, 0, 1)
            },
            "v-shape" = {
              p_val <- if (is.list(thresh) && !is.null(thresh$p)) thresh$p else if (is.numeric(thresh) && length(thresh) == 1) thresh else stop(paste("Missing or invalid 'p' threshold for v-shape preference function (criterion", k, ")"))
              ifelse(d <= 0, 0, d / p_val)
            },
            "level" = {
              if (is.list(thresh) && !is.null(thresh$p) && !is.null(thresh$q)) {
                p_val <- thresh$p
                q_val <- thresh$q
              } else if (is.numeric(thresh) && length(thresh) == 2) {
                p_val <- thresh[1]
                q_val <- thresh[2]
              } else {
                stop(paste("Missing 'p' or 'q' thr for level preference function (criterion", k, ")"))
              }
              ifelse(d <= q_val, 0, ifelse(d > p_val, 1, 0.5))
            },
            stop(paste("Invalid preference function:", pref_type, "(criterion", k, ")"))
          )
        }
      }
    }
  }

  # Weighted Preference Indices
  weighted_preference_indices <- apply(preference_indices, c(1, 2), function(x) sum(x * weights))

  # Calculate Positive and Negative Outranking Flows
  phi_plus <- rowSums(weighted_preference_indices) / (n - 1)
  phi_minus <- colSums(weighted_preference_indices) / (n - 1)

  # Calculate Net Outranking Flow (PROMETHEE II)
  phi_net <- phi_plus - phi_minus

  # Central Tendency and Variability for PROMETHEE VI
  if (is.null(p)) {
    warning("p parameter is NULL. Default value is assigned based on the range of criterion values.")
    p <- max(apply(dmatrix, 2, function(x) diff(range(x)))) # Maximum of the range of each criterion
  }
  if (is.null(q)) {
    warning("q parameter is NULL. Default value is assigned based on the standard deviation of criterion values.")
    q <- mean(apply(dmatrix, 2, sd)) / 2 # Half of the average standard deviation of the criteria
  }
  if (!is.numeric(p) || !is.numeric(q) || length(p) != 1 || length(q) != 1 || p <= q) {
    stop("'p' and 'q' must be numeric, single values, and 'p' > 'q'.")
  }

  central_tendency <- phi_net # Net flow as central tendency

  # Calculate Variability
  variability <- switch(varmethod,
    "abs_sum" = abs(phi_plus + phi_minus), # Absolute sum (recommended)
    "range" = phi_plus + phi_minus,       # Flow range (alternative)
    stop("Invalid variability_method. Must be 'abs_sum' or 'range'.")
  )

  # Calculate Lambda
  lambda <- numeric(n)
  for (i in 1:n) {
    if (variability[i] <= q) {
      lambda[i] <- 1
    } else if (variability[i] > q && variability[i] < p) {
      lambda[i] <- (p - variability[i]) / (p - q)
    } else {
      lambda[i] <- 0
    }
  }

  # Calculate the final value Xi
  xi <- central_tendency + lambda * variability

  # Ranking
  final_ranking <- rank(-xi, ties.method = tiesmethod)
  names(final_ranking) <- rownames(dmatrix)
  
  results <- list(
    Phi_plus = phi_plus,
    Phi_minus = phi_minus,
    Phi_net = phi_net,
    Central_tendency = central_tendency,
    Variability = variability,
    Lambda = lambda,
    Xi = xi,
    rank = final_ranking
  )
  return(results)
}