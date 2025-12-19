ftopsis <- function(fuzzydmatrix, bcvec, fuzzyweights, tiesmethod = "average") {

  # Input Checks 
  # fuzzydmatrix must be a matrix where each fuzzy number (l, m, u) occupies 3 columns.
  if (!is.matrix(fuzzydmatrix) || ncol(fuzzydmatrix) %% 3 != 0) {
    stop("Fuzzy decision matrix (fuzzydmatrix) must be a matrix, and its number of columns must be a multiple of 3 (for l, m, u format).")
  }

  m <- ncol(fuzzydmatrix) / 3

  if (length(bcvec) != m) {
    stop("The benefit/cost vector (bcvec) must have the same number of elements as the criteria.")
  }

  if (!all(bcvec %in% c(1, -1))) {
    stop("The benefit/cost vector (bcvec) must contain 1 (benefit) or -1 (cost) values.")
  }

  # fuzzyweights must be a matrix with 'm' rows and 3 columns (l, m, u).
  if (!is.matrix(fuzzyweights) || ncol(fuzzyweights) != 3 || nrow(fuzzyweights) != m) {
    stop("Weights (fuzzyweights) must be a matrix with 'm' rows and 3 columns (l, m, u).")
  }

  # Fuzzy Normalization 
  # Normalization for triangular fuzzy numbers.
  # For benefit criteria: X_ij / u_j_max
  # For cost criteria: (l_j_min / u_ij, l_j_min / m_ij, l_j_min / l_ij)
  
  nfmatrix <- matrix(NA, nrow = nrow(fuzzydmatrix), ncol = ncol(fuzzydmatrix))
  
  for (j in 1:m) {
    # Column indices for the current criterion's fuzzy number (l, m, u)
    cidxl <- (j - 1) * 3 + 1
    cidxm <- (j - 1) * 3 + 2
    cidxu <- (j - 1) * 3 + 3
    
    current_criterion_fuzzy_vals <- fuzzydmatrix[, c(cidxl, cidxm, cidxu)]

    if (bcvec[j] == 1) { # Benefit criterion
      max_u_val <- max(current_criterion_fuzzy_vals[, 3]) # Max 'u' value across alternatives for this criterion
      # Divide each component (l, m, u) by max_u_val
      nfmatrix[, cidxl] <- current_criterion_fuzzy_vals[, 1] / max_u_val
      nfmatrix[, cidxm] <- current_criterion_fuzzy_vals[, 2] / max_u_val
      nfmatrix[, cidxu] <- current_criterion_fuzzy_vals[, 3] / max_u_val     
    } else { # Cost criterion
      min_l_val <- min(current_criterion_fuzzy_vals[, 1]) # Min 'l' value across alternatives for this criterion
      # Normalize for cost criterion: (min_l / u, min_l / m, min_l / l)
      nfmatrix[, cidxl] <- min_l_val / current_criterion_fuzzy_vals[, 3] # New l = min_l / original u
      nfmatrix[, cidxm] <- min_l_val / current_criterion_fuzzy_vals[, 2] # New m = min_l / original m
      nfmatrix[, cidxu] <- min_l_val / current_criterion_fuzzy_vals[, 1] # New u = min_l / original l
    }
  }

  # Create Weighted Fuzzy Normalized Decision Matrix 
  # Fuzzy multiplication: (l1*l2, m1*m2, u1*u2)
  wfmatrix <- matrix(NA, nrow = nrow(fuzzydmatrix), ncol = ncol(fuzzydmatrix))

  for (j in 1:m) {
    cidxl <- (j - 1) * 3 + 1
    cidxm <- (j - 1) * 3 + 2
    cidxu <- (j - 1) * 3 + 3

    # Multiply normalized fuzzy values by corresponding fuzzy weights
    wfmatrix[, cidxl] <- nfmatrix[, cidxl] * fuzzyweights[j, 1]
    wfmatrix[, cidxm] <- nfmatrix[, cidxm] * fuzzyweights[j, 2]
    wfmatrix[, cidxu] <- nfmatrix[, cidxu] * fuzzyweights[j, 3]
  }

  # Determine Fuzzy Positive Ideal Solution (FPIS) and Fuzzy Negative Ideal Solution (FNIS) 
  FPIS <- numeric(ncol(fuzzydmatrix)) # Initialize as numeric vector
  FNIS <- numeric(ncol(fuzzydmatrix)) # Initialize as numeric vector

  for (j in 1:m) {
    cidxl <- (j - 1) * 3 + 1
    cidxm <- (j - 1) * 3 + 2
    cidxu <- (j - 1) * 3 + 3
    
    # FPIS: Max 'l', 'm', 'u' for each criterion column
    FPIS[cidxl] <- max(wfmatrix[, cidxl])
    FPIS[cidxm] <- max(wfmatrix[, cidxm])
    FPIS[cidxu] <- max(wfmatrix[, cidxu])
    
    # FNIS: Min 'l', 'm', 'u' for each criterion column
    FNIS[cidxl] <- min(wfmatrix[, cidxl])
    FNIS[cidxm] <- min(wfmatrix[, cidxm])
    FNIS[cidxu] <- min(wfmatrix[, cidxu])
  }

  # Calculate Distances from FPIS and FNIS (using Vertex Method) 
  # Distance d(A, B) = sqrt( (l_A-l_B)^2 + (m_A-m_B)^2 + (u_A-u_B)^2 )
  
  d_plus_fuzzy <- numeric(nrow(fuzzydmatrix))  # Distance to FPIS for each alternative
  d_minus_fuzzy <- numeric(nrow(fuzzydmatrix)) # Distance to FNIS for each alternative

  for (i in 1:nrow(fuzzydmatrix)) {
    sum_sq_dist_plus <- 0 # Sum of squared differences for distance to FPIS
    sum_sq_dist_minus <- 0 # Sum of squared differences for distance to FNIS
    
    for (j in 1:m) {
      cidxl <- (j - 1) * 3 + 1
      cidxm <- (j - 1) * 3 + 2
      cidxu <- (j - 1) * 3 + 3
      
      # Get the fuzzy value for the current alternative and criterion
      alt_l <- wfmatrix[i, cidxl]
      alt_m <- wfmatrix[i, cidxm]
      alt_u <- wfmatrix[i, cidxu]
      
      # Add squared differences for FPIS distance
      sum_sq_dist_plus <- sum_sq_dist_plus +
          (alt_l - FPIS[cidxl])^2 +
          (alt_m - FPIS[cidxm])^2 +
          (alt_u - FPIS[cidxu])^2
      
      # Add squared differences for FNIS distance
      sum_sq_dist_minus <- sum_sq_dist_minus +
          (alt_l - FNIS[cidxl])^2 +
          (alt_m - FNIS[cidxm])^2 +
          (alt_u - FNIS[cidxu])^2
    }

    d_plus_fuzzy[i] <- sqrt(sum_sq_dist_plus)
    d_minus_fuzzy[i] <- sqrt(sum_sq_dist_minus)
  }

  # Calculate Closeness Coefficient (CCi) 
  # CCi = D_minus / (D_minus + D_plus)
  # Handle potential division by zero if D_minus + D_plus is very small or zero

  denominator <- d_minus_fuzzy + d_plus_fuzzy
  closeness_coefficient_fuzzy <- ifelse(denominator == 0, 0, d_minus_fuzzy / denominator)

  # Higher CCi value indicates a better alternative.
  ranking <- rank(-closeness_coefficient_fuzzy, ties.method = tiesmethod)

  results <- list(
    Criteria = m,
    fuzzyweights = fuzzyweights,
    FPIS = FPIS,
    FNIS = FNIS,
    distance_to_FPIS = d_plus_fuzzy,
    distance_to_FNIS = d_minus_fuzzy,
    closeness_coefficient_fuzzy = closeness_coefficient_fuzzy,
    rank = ranking
  )

  return(results)
}
