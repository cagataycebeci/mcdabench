# Multi-criteria Evaluation with Gradually-weighted Aggregation of Numerical Decision Matrices
megan2 <- function(dmatrix, bcvec, weights, weimethod="random", normethod=NULL, thr=0, tht="p25", 
  tiesmethod="average", ss=0.05, niters=10, rp=c(0.1, 0.5, 0.1)){

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  if(missing(weights)){   
     weights <- rep(1/m, m)
  }
  
  # Gradually modify the weights
  if(weimethod=="gradual"){
     weimat <- gengradwei(weights, rp=rp)
  } else if(weimethod=="random"){
     weimat <- genrandwei(weights, ss=ss, niters=niters)
  } else {
    stop("Invalid weight modification method. Set the weimethod argument to 'gradual' or 'random'.")
  }
  
  rankmat <- matrix(NA, nrow=nrow(weimat), ncol=n, dimnames=list(rownames(weimat), rownames(dmatrix)))
  
  for(i in 1:nrow(rankmat)){
   rankmat[i,] <- megan(dmatrix, bcvec=bcvec, weights=weimat[i,], normethod=normethod, thr=thr, tht=tht)$rank
  }

  #ranking <- rankaggregate(rankmat, tiesmethod=tiesmethod)$preference_ranking["BORDACNT",]
  # Borda Count ranking
 
  borda_scores <- numeric(ncol(rankmat)) 
  names(borda_scores) <- colnames(rankmat)  
    
  for (i in 1:nrow(rankmat)) {  
      borda_scores <- borda_scores + (nrow(rankmat) - rankmat[i, ])
  }
    
  borda_ranking <- rank(-borda_scores, ties.method = tiesmethod) 
     
  results <- list(weimat=weimat, 
     rankmat=rankmat, 
     rank=borda_ranking
  )
  
  return(results)
}
