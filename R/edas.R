# Evaluation-based on Distance from Average Solution (EDAS) 
edas <- function(dmatrix, bcvec, weights, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop(" The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  pda <- matrix(0, nrow = n, ncol = m)
  nda <- matrix(0, nrow = n, ncol = m)
  
  avec <- colMeans(dmatrix)

  for (j in 1:m) {
    if (bcvec[j] == -1) {
      pda[, j] <- (avec[j] - dmatrix[, j]) / avec[j]
      nda[, j] <- (dmatrix[, j] - avec[j]) / avec[j]
    } else {
      pda[, j] <- (dmatrix[, j] - avec[j]) / avec[j]
      nda[, j] <- (avec[j] - dmatrix[, j]) / avec[j]
    }
  }
  
  pda[pda < 0] <- 0
  nda[nda < 0] <- 0
  
  sp <- rowSums(weights * pda)
  sn <- rowSums(weights * nda)
  
  nsp <- sp / max(sp)
  nsn <- 1 - sn / max(sn)
  
  score <- (nsp + nsn) / 2
  
  ranking <- rank(-score, ties.method=tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  result <- list(
    avec = avec, 
    pda = pda, 
    nda = nda, 
    sp = sp, 
    sn = sn, 
    nsp = nsp, 
    nsn = nsn, 
    score = score, 
    rank=ranking
  )
  
  return(result)
}

