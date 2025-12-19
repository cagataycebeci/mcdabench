smart <- function(dmatrix, bcvec=NULL, weights, normethod=NULL, valfuncs = NULL, tiesmethod="average") {
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
  
  if(is.null(normethod)){
    nmatrix <- dmatrix
  }else{
    # Normalize the decision matrix
    nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type=normethod)
  }
  
  w <- weights / sum(weights) # calcnormal weights

  # Define value functions (default to linear)
  if (is.null(valfuncs)) {
    valfuncs <- lapply(1:m, function(j) {
      min_val <- min(nmatrix[, j])
      max_val <- max(nmatrix[, j])
      if (min_val == max_val) {
        return(function(x) 0.5) 
      } else {
        return(function(x) (x - min_val) / (max_val - min_val))
      }
    })
  } else if (!is.list(valfuncs) || length(valfuncs) != m) {
    stop("Value functions (valfuncs) must be a list of functions with length equal to the number of columns in dmatrix.")
  }

  # Calculate the value matrix
  vmatrix <- matrix(0, nrow = n, ncol = m)
  for (j in 1:m) {
    vmatrix[, j] <- sapply(nmatrix[, j], valfuncs[[j]])
  }

  # Calculate weighted value matrix
  wvmatrix <- vmatrix * w

  # Calculate overall scores
  overall_scores <- rowSums(wvmatrix)

  ranking <- rank(-overall_scores, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
 
  results <- list(
    Value_Matrix = vmatrix,
    Weighted_Value_Matrix = wvmatrix,
    Overall_Scores = overall_scores,
    rank = ranking
  )
  
  return(results)
}