mavt <- function(dmatrix, bcvec, weights, valfuncs = NULL, normvals = TRUE, ss = 1, tiesmethod="average") {

  # Define Marginal Utility Functions
  # Exponential
  vfexp <- function(x) {
    y <- (exp(x^2) - 1) / 1.72
    return(y)
  }

  # Log10-based
  vflog <- function(x) {
    y <- log10(9 * x + 1)
    return(y)
  }

  # Natural Logarithm-based
  vfln <- function(x) {
    y <- log((exp(1) - 1) * x + 1)
    return(y)
  }

  # Quadratic
  vfquad <- function(x) {
    y <- (2 * x - 1)^2
    return(y)
  }
  
  # Step with a step size
  vfstep <- function(x, ss) {
    y <- ceiling(ss * x) / ss
    return(y)
  }
  
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("dmatrix must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec (benefit/cost vector) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  # Default to linear value functions
  if (is.null(valfuncs)) {
    valfuncs <- rep("linear", m)
  } else if (length(valfuncs) != m && length(valfuncs) != 1) {
    stop("The length of the value functions vector must be 1 or equal to the number of criteria.")
  } else if (length(valfuncs) == 1) {
    valfuncs <- rep(valfuncs, m)
  }

  value_matrix <- matrix(0, nrow = n, ncol = m)
  colnames(value_matrix) <- colnames(dmatrix)
  rownames(value_matrix) <- rownames(dmatrix)

  for (j in 1:m) {
    values <- dmatrix[, j]
    if (valfuncs[j] == "linear") {
      vj <- values
    } else if (valfuncs[j] == "exp") {
      vj <- vfexp(values)
    } else if (valfuncs[j] == "log") {
      vj <- vflog(values)
    } else if (valfuncs[j] == "ln") {
      vj <- vfln(values)
    } else if (valfuncs[j] == "quad") {
      vj <- vfquad(values)
    } else if (valfuncs[j] == "step") {
      vj <- vfstep(values, ss)
    } else {
      stop(paste("Unrecognized value function:", valfuncs[j]))
    }    
    value_matrix[, j] <- vj
  }

  # Optionally calcnormal the value matrix to a 0-1 range
  if (normvals) {
    value_matrix <- apply(value_matrix, 2, function(x) {
      min_x <- min(x)
      max_x <- max(x)
      if (max_x == min_x) {
        return(rep(0.5, length(x))) 
      } else {
        return((x - min_x) / (max_x - min_x))
      }
    })
  }
  
  # Calculate the weighted average value scores
  overall_scores <- rowSums(value_matrix * weights)
  names(overall_scores) <- rownames(dmatrix)
  
  ranking <- rank(-overall_scores, ties.method = tiesmethod)
  names(ranking) <- names(overall_scores)  
  
  results <- list(
    Value_Matrix = value_matrix,
    Weighted_Sum_Scores = overall_scores,
    rank = ranking
  )
  return(results)
}
