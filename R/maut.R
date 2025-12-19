# Multi-Attribute Utility Theory (MAUT) evaluation. 
maut <- function(dmatrix, bcvec, weights, utilfuncs = NULL, normutil = TRUE, ss=1, tiesmethod = "average") {

  # Define Marginal Utility Functions
  # Exponential
  ufexp <- function(x) {
    y <- (exp(x^2) - 1) / 1.72
    return(y)
  }

  # Log10-based
  uflog <- function(x) {
    y <- log10(9 * x + 1)
    return(y)
  }

  # Natural Logarithm-based
  ufln <- function(x) {
    y <- log((exp(1) - 1) * x + 1)
    return(y)
  }

  # Quadratic
  ufquad <- function(x) {
    y <- (2 * x - 1)^2
    return(y)
  }
  
  # Step with a step size
  ufstep <- function(x, ss) {
    y <- ceiling(ss * x) / ss
    return(y)
  }

  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (for cost) or 1 (for benefit) values with a length equal to the number of criteria.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  # Set default utility functions to "linear" if not provided
  if (is.null(utilfuncs)) {
    utilfuncs <- rep("linear", m)
  } else if (length(utilfuncs) != m && length(utilfuncs) != 1) {
    stop("The length of the utility functions vector must be 1 or equal to the number of criteria.")
  } else if (length(utilfuncs) == 1) {
    utilfuncs <- rep(utilfuncs, m)
  }
  
  # Create a utility matrix of the same size as the decision matrix
  utility_matrix <- matrix(0, nrow = n, ncol = m)
  colnames(utility_matrix) <- colnames(dmatrix)
  rownames(utility_matrix) <- rownames(dmatrix)
  
  # Apply min-max normalization and then the chosen utility function for each criterion
  for (j in 1:m) {
    values <- dmatrix[, j]
    
    # Perform min-max normalization based on criterion type:
    if (bcvec[j] == 1) {  
      normval <- (values - min(values)) / (max(values) - min(values))
    } else {  # Cost criterion (bcvec = -1)
      normval <- (max(values) - values) / (max(values) - min(values))
    }
    
    # Apply the specified utility function transformation:
    if (utilfuncs[j] == "linear") {
      uj <- normval
    } else if (utilfuncs[j] == "exp") {
      uj <- ufexp(normval)
    } else if (utilfuncs[j] == "log") {
      uj <- uflog(normval)
    } else if (utilfuncs[j] == "ln") {
      uj <- ufln(normval)
    } else if (utilfuncs[j] == "quad") {
      uj <- ufquad(normval)
    } else if (utilfuncs[j] == "step") {
       uj <- ufstep(normval, ss)
    } else {
      stop(paste("Unrecognized utility function:", utilfuncs[j]))
    }
    
    utility_matrix[, j] <- uj
  }
  
  # Optionally, re-normalze the utility matrix to the [0,1] range column-wise
  if (normutil) {
    utility_matrix <- apply(utility_matrix, 2, function(x) {
      if (diff(range(x)) != 0) {
        (x - min(x)) / (max(x) - min(x))
      } else {
        rep(0.5, length(x))  # If values are constant, assign a mid-range value
      }
    })
  }
  
  # Calculate the weighted average utility scores for each alternative
  overall_scores <- rowSums(utility_matrix * weights)
  ranking <- rank(-overall_scores, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  results <- list(
    Utility_Matrix = utility_matrix,
    Weighted_Sum_Scores = overall_scores,
    rank = ranking
  )
  
  return(results)
}
