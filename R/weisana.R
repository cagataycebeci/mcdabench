weisana <- function(dmatrix, bcvec, weights, 
           weimethod = "gradual", weipars = NULL,
           mcdamethod = topsis, methodpars, sensplot=TRUE) { 
  
  if(is.null(weipars)){
     weipars <- list(
       rp = seq(0.1, 0.5, 0.1),
       ss = 0.01,
       niters = 10
     )
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  if(missing(weights)){   
       weights <- rep(1/m, m)
  }
  
  if(missing(methodpars)){   
     methodpars <- list()
     #stop("Define the parameters of method with methodpars")
  }
  
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
 
  # Generate modified weights
  if(weimethod %in% c("g","grad", "gradual")){
    weimat <- gengradwei(weights, rp = weipars$rp)
  } else if (weimethod %in% c("r","rand", "random")) {
    weimat <- genrandwei(weights, ss = weipars$ss, niters=weipars$niters)
  } else {
    stop("Invalid weight modification method. Set the weimethod argument to 'gradual' or 'random'.")
  }
    
  rankmat <- matrix(NA, nrow=nrow(weimat), ncol=n, dimnames=list(rownames(weimat), rownames(dmatrix)))
  colnames(rankmat) <- rownames(dmatrix)
  
  starttime <- Sys.time()
  for(i in 1:nrow(rankmat)){
     args <- c(list(dmatrix = dmatrix, bcvec = bcvec, weights = weimat[i,]), methodpars)
     resiter <- do.call(mcdamethod, args)
     rankmat[i,] <- resiter$rank
  }
  endtime <- Sys.time()
  comptime <- round(difftime(endtime, starttime, units = "secs"), 2)
  
  n <- nrow(rankmat)
  uniqueranks <- as.matrix(unique(rankmat))
  resunique <- data.frame(
    Pattern = apply(uniqueranks, 1, paste, collapse = ","), 
    Count = NA, Percent = NA, stringsAsFactors = FALSE)
  rownames(resunique) <- paste0("Ranking_", 1: nrow(resunique))
  for (i in 1:nrow(uniqueranks)) {
    cuniqpat <- uniqueranks[i, ]
    matches <- apply(rankmat, 1, function(x) {
        all(x == cuniqpat) 
    })
    count <- sum(matches)
    resunique$Count[i] <- count
    resunique$Percent[i] <- round((count / n) * 100, 2)
  } 
  
  if(sensplot) sensplot(resunique)
  
  results <- list(
     computing_time=comptime, 
     weights_matrix=weimat, 
     ranking_matrix=rankmat, 
     sensitivity_table=resunique
  )

  return(results)
}
