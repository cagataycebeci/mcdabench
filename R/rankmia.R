# Modified Index of Agreement
rankmia <- function(rankmat, j = 1, na.rm = TRUE) {
  n <- nrow(rankmat)
  rnames <- rownames(rankmat)
  ranklist <- lapply(1:n, function(i) rankmat[i, ])
  names(ranklist) <- rnames

  midmatrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))

  midx <- function(sim, obs, j = 1, na.rm = TRUE) {
    if (na.rm) {
      ind <- !(is.na(sim) | is.na(obs))
      sim <- sim[ind]
      obs <- obs[ind]
    }
    nval <- length(obs)
    if (nval == 0) return(NA)

    numerator <- sum(abs(obs - sim)^j)
    denominator <- sum((abs(sim - mean(obs))^j) + (abs(obs - mean(obs))^j))

    if (denominator == 0) {
      return(1) 
    } else {
      return(1 - (numerator / denominator))
    }
  }

  for (i in 1:n) {
    for (k in 1:n) {
      midval <- midx(ranklist[[i]], ranklist[[k]], j = j, na.rm = na.rm)
      midmatrix[i, k] <- midval
    }
  }
   
  displaymat <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))
  for (i in 1:n) {
      for (j in 1:n) {
        if (i < j) {
          displaymat[i, j] <- format(round(midmatrix[i, j], 2), nsmall = 2)
        } else if (i > j) {
          displaymat[i, j] <- format(round(midmatrix[i, j], 2), nsmall = 2)
        } else {
          displaymat[i, j] <- "" 
        }
      }
    }
  displaymat <- as.data.frame(displaymat)
  
  result <- list(midmatrix=midmatrix, displaymat=displaymat)
  return(result)
}