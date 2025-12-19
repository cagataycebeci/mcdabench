# Measurement of Alternatives and Ranking According to Compromise Solution (MARCOS)
marcos <- function(dmatrix, bcvec, weights, normethod = "ratio", tiesmethod="average") {
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
  
  #Extended decision matrix
  ematrix <- matrix(0, n + 2, m)
  ematrix[1:n, ] <- dmatrix
  
  colmax <- apply(dmatrix, 2, max)
  colmin <- apply(dmatrix, 2, min)
   
  for (i in 1:m) {
    if (bcvec[i] == 1) {
      ematrix[n + 1, i] <- colmax[i]
      ematrix[n + 2, i] <- colmin[i]
    } else {
      ematrix[n + 1, i] <- colmin[i]
      ematrix[n + 2, i] <- colmax[i]
    }
   }
  
   nexmatrix <- calcnormal(ematrix, bcvec=bcvec, type=normethod)
   wmatrix <- sweep(nexmatrix, 2, weights, "*")

   S <- rowSums(wmatrix)
   # Check for division by zero
   if (S[n + 2] == 0) S[n + 2] <- 1e-6 
   if (S[n + 1] == 0) S[n + 1] <- 1e-6 

   # Utility calculations
   k_neg <- S[1:n] / S[n + 2]
   k_pos <- S[1:n] / S[n + 1]
   
   f_k_pos <- k_neg / (k_pos + k_neg)
   f_k_neg <- k_pos / (k_pos + k_neg)
   f_k <- (k_pos + k_neg) / (1 + (1 - f_k_pos) / f_k_pos + (1 - f_k_neg) / f_k_neg)
  
   ranking <- rank(-f_k, ties.method=tiesmethod)
   names(ranking) <- rownames(dmatrix)
   
   result <- list(ematrix = ematrix, nexmatrix = nexmatrix, 
     wmatrix = wmatrix, S = S, k_neg = k_neg, k_pos = k_pos, 
     f_kpos = f_k_pos, f_kneg = f_k_neg, f_k = f_k, rank=ranking)
  
   return(result)
}
