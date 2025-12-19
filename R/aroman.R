# Alternative Ranking Order Method Accounting for Two Step Normalization (AROMAN)
aroman <- function(dmatrix, bcvec, weights, beta=0.5, lambda=0.5, tiesmethod="average") {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if(beta < 0 || beta> 1) 
    stop("Beta must be between 0 and 1.")
 
  if(lambda < 0 || lambda > 1) 
    stop("Lambda must be between 0 and 1.")
  
  m <- nrow(dmatrix) 
  n <- ncol(dmatrix) 

  nmatrix1 <- calcnormal(dmatrix, bcvec, type="maxmin")
  nmatrix2 <- calcnormal(dmatrix, bcvec, type="vector")

  Tnorm <- (beta * nmatrix1 + (1 - beta) * nmatrix2) / 2

  if(missing(weights)){
    W <- rep(1/n, n)
  } else {
    if(length(weights) != n) {
      stop("Check the lenght of weights vector.")
    }
    W <- weights
  }

  if (!is.numeric(W) || length(W) != ncol(dmatrix)) {
      stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }

  That <- Tnorm * matrix(W, nrow = m, ncol = n, byrow = TRUE)

  L <- numeric(m) # For cost criteria
  A <- numeric(m) # For benefit criteria

  if(length(bcvec) != n) {
    stop("Check the lenght of bcvec!")
  }

  for (i in 1:m) {
    L[i] <- sum(That[i, which(bcvec == -1)])
    A[i] <- sum(That[i, which(bcvec == 1)])
  }

  R <- (L^lambda + A^(1 - lambda))^(1 / (1 - lambda))

  #ranking <- order(R, decreasing = TRUE)
  ranking <- rank(-R, ties.method="average")
 
  results <- list(
    L = L,
    A = A,
    R = R,
    rank = ranking
  )
  return(results)
}
