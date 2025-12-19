# Root Assessment Method (RAM)
ram <- function(dmatrix, bcvec, weights, normethod = "sum", tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
 
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Normalized the decision matrix
  nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type=normethod)
  
  # Perform element-wise multiplication with weights
  wmatrix <- sweep(nmatrix, 2, weights, "*")
  
  # Mask for benefit-cost vector
  benefit_mask <- which(bcvec == 1)
  cost_mask <- which(bcvec == -1)
  
  # Fix rowSums issue by ensuring proper matrix structure
  Spi <- rowSums(as.matrix(wmatrix[, benefit_mask, drop = FALSE]))
  Smi <- rowSums(as.matrix(wmatrix[, cost_mask, drop = FALSE]))
  
  # Calculate the root assessment index
  ri <- (2 + Spi) ^ (1 / (2 + Smi))
  
  ranking <- rank(-ri, ties.method=tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  result <- list(
    nmatrix = nmatrix,
    weighted_matrix = wmatrix,
    Spi = Spi,
    Smi = Smi,
    ri = ri,
    rank = ranking
  )
  
  return(result)
}
