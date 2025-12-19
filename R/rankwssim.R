rankwssim <- function(rankmat) {
  # WS Coefficient calculation
  wsrc <- function(rnk1, rnk2) {
    n <- length(rnk1)
    total <- 0
    for (i in 1:n) {
      R_y_i <- rnk2[i] 
      R_x_i <- rnk1[i] 
      
      nominator <- abs(R_x_i - R_y_i)
      denominator <- max(abs(1 - R_x_i), abs(n - R_x_i)) + 1e-10  

      total <- total + (2^(-R_x_i) * (nominator / denominator))
    }
    score <- 1 - total
    return(max(0, min(score, 1)))
  }

  n <- nrow(rankmat) 
  rnames <- rownames(rankmat)
  ranklist <- lapply(1:n, function(i) rankmat[i, ])
  names(ranklist) <- rnames

  wsmatrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))

  for (i in 1:n) {
    for (j in 1:n) {
      wsmatrix[i, j] <- wsrc(ranklist[[i]], ranklist[[j]])
    }
  }

  displaymat <- matrix("", nrow = n, ncol = n, dimnames = list(rnames, rnames))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        displaymat[i, j] <- format(round(wsmatrix[i, j], 2), nsmall = 2)
      } else if (i > j) {
        displaymat[i, j] <- format(round(wsmatrix[i, j], 2), nsmall = 2)
      } else {
        displaymat[i, j] <- "" 
      }
    }
  }
  displaymat <- as.data.frame(displaymat)
  
  return(list(
    wsmatrix = wsmatrix, 
    displaymat = displaymat 
  ))
}
