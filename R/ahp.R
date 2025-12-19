ahp <- function(pmatrix, method="maxeig", weights=NULL, tiesmethod="average") {
  # Predefined Random Index values from Saaty's AHP model
  ri <- c(0, 0, 0, 0.58, 0.9, 1.12, 1.24, 1.32, 1.41,
          1.45, 1.49, 1.51, 1.48, 1.56, 1.57, 1.59)
  
  k <- length(ri)

  pmatrix <- as.matrix(pmatrix)
  n <- ncol(pmatrix)

  if (nrow(pmatrix) != ncol(pmatrix)) {
    stop("Error: The comparison matrix must be square (same number of rows & columns).")
  }

  if (any(pmatrix <= 0)) {
    stop("Error: All elements in the comparison matrix must be positive.")
  }
  if (!all(diag(pmatrix) == 1)) {
    stop("Error: All diagonal elements must be 1.")
  }
  if (!isTRUE(all.equal(pmatrix, 1/t(pmatrix), tolerance = 1e-9))) {
    stop("Error: The comparison matrix must be reciprocal (a_ij must be 1/a_ji).")
  }

  # Extend Random Index values for large matrices using linear interpolation
  extri <- if (n > k) {
    ri_extended <- ri 
    # Loop from the first value beyond the predefined table up to 'n'
    for (i in (k + 1):n) {
        # Calculate the slope using the last two known/calculated RI values
        # The indices for lookup are 'i-1' and 'i-2'
        slope <- (ri_extended[i-1] - ri_extended[i-2]) / ((i-1) - (i-2))       
        # Extrapolate the next RI value
        ri_extended[i] <- ri_extended[i-1] + slope * (i - (i-1))
    }
    ri_extended[1:n]
  } else {
    ri
  }
 
  # If user provides weights, use them directly
  if (!is.null(weights)) {
    weights <- weights / sum(weights)  # Normalize weights
  } else {
    # Compute weights based on method
    if (method %in% c("m", "mean", "arit", "arithmetic")) {
      weights <- apply(pmatrix, 1, function(x) mean(x / colSums(pmatrix)))
    } else if (method %in% c("g", "geo", "geom", "geometric")) {
      weights <- apply(pmatrix, 1, function(x) prod(x)^(1/n))
      weights <- weights / sum(weights)
    } else if (method %in% c("me", "maxeig")) {
      eig <- eigen(pmatrix)
      preigvec <- Re(eig$vectors[, which.max(Re(eig$values))])
      weights <- preigvec / sum(preigvec)
    }
  }
  
  # Compute maximum eigenvalue & consistency index
  lambdamax <- mean(rowSums(pmatrix * weights) / weights)
  idxcons <- (lambdamax - n) / (n - 1)

  # Calculate Consistency Ratio safely
  safe_extri <- extri[min(n, length(extri))]
  rc <- ifelse(safe_extri <= 0, NA, idxcons / safe_extri)

  ranking <- rank(-weights, ties.method = tiesmethod)

  results <- list(weights = weights, rc = rc, rank = ranking)
  
  return(results)
}
