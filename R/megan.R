# Multi-criteria Evaluation with MEGAN
megan <- function(dmatrix, bcvec, weights="gini", normethod=NULL, thr=NULL, tht="p25", tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("dmatrix must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("bcvec (benefit/cost vector) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  
  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  # Normalize the decision matrix
  if(is.null(normethod)){
    nmatrix <- dmatrix
  } else {
    nmatrix <- calcnormal(dmatrix, bcvec=bcvec, type=normethod)
  }
  
  # Calculate the weights of criteria
  if (is.character(weights)) {
      w <- calcweights(nmatrix, bcvec, type = weights) 
  } else if (is.numeric(weights) && length(weights) > 1) {
      w <- weights
  } 
  if (!is.numeric(w) || length(w) != m) {
     stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }
    
  # Perform element-wise multiplication with weights
  wmatrix <- sweep(nmatrix, 2, w, "*") 

  # Convert the costs to benefits
  amatrix <- sweep(wmatrix, 2, bcvec, "*") 
  
  # Find the minimum values of criteria
  minvec <- abs(apply(amatrix, 2, min))
  
  # Calculate distances to minimum value
  distmatrix <- sweep(amatrix, 2, minvec, "-")
  
  # calcnormal distance matrix
  ndistmatrix <- apply(distmatrix, 2, FUN=function(x) (x-min(x)) / (max(x)-min(x)))

  # Calculate row sums of calcnormald distances
  distrowsums <- apply(ndistmatrix, 1, sum) 

  # Calculate superiority scores of alternatives (as percentage over the minimum)
  if(min(distrowsums) == 0) {
      superiority <- round((distrowsums - min(distrowsums)) * 100)
  } else {
    superiority <- round((distrowsums-min(distrowsums))/min(distrowsums) * 100)
  }
  # Determine threshold for ranking if not supported
  if(is.null(thr)){
    if(tht == "sdev"){
       thr <- sd(diff(sort(superiority))) / 2
    } else if(tht == "p5"){
       thr <- quantile(diff(sort(superiority)), 0.05)
    } else if(tht == "p25"){
       thr <- quantile(diff(sort(superiority)), 0.25)
    } else {
      stop("Invalid threshold type for ranking. It should be  one of 'sd', 'p5', 'p25'.")
    }
  } else {
    # Check the validity of threshold
    if(!is.numeric(thr) || (thr < 0 || thr > 100)) {
      stop("Threshold value should be a positive number between 0 and 100. Please use 'p5', 'p25' or 'sdev' for computing it internally.")
    }
    tht <- "user"
  }
  
  # Find the superiority ranks using threshold
  sidx <- order(-superiority)  
  ssuperiority <- superiority[sidx]  
  sranks <- numeric(length(superiority))  
  sranks[1] <- 1  
  
  for (i in 2:length(superiority)) {  
    if (abs(ssuperiority[i] - ssuperiority[i - 1]) > thr) {  
      sranks[i] <- sranks[i - 1] + 1  
    } else {  
      sranks[i] <- sranks[i - 1]  
    }  
  }  
  
  sranks_final <- numeric(length(superiority))  
  sranks_final[sidx] <- sranks  
  
  rank <- rank(sranks_final, ties.method=tiesmethod)  
  names(rank) <- rownames(dmatrix)  
  
  result <- list(dmatrix=dmatrix, wmatrix=wmatrix, minvec=minvec,  
     distmatrix=distmatrix, ndistmatrix=ndistmatrix, distrowsums=distrowsums,  
     superiority=superiority, threshmethod=tht, thresholdvalue=thr, rank=rank)  
  
  return(result)
}
