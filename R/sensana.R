# Sensitivity and Stability Analyses for MCDA
sensana <- function(rankmat){
  n <- nrow(rankmat)
  m <- ncol(rankmat)

  # Standard deviation
  sdev <- apply(rankmat, 2, sd)

  # Standard deviation two stages
  minrank <- min(rankmat)
  maxrank <- max(rankmat)
  nrankmat <- (rankmat - minrank) / (maxrank - minrank)
  sdev2 <- apply(nrankmat, 2, function(x) mean(abs(x - mean(x))))
 
  # Coefficient of Rank Variation (CVR)
  crv <- apply(rankmat, 2, function(x) {sd(x) / mean(x)})

  # Rank Stability Index (aka Stability scores)
  stabscores <- numeric(m)

  for (j in 1:m) {
    ranks <- rankmat[,j]
    uniqueranks <- length(unique(ranks))
    
    # Stability score: 1 - (percentage of changing ranks)
    stabscores[j] <- 1 - ((uniqueranks - 1) / n)
  }

  # Sequential Rank Shift Index (SRSI)
  srsi <- numeric(m)  
  for (j in 1:m) {
    ranks <- rankmat[,j]  
    if(length(unique(ranks)) == 1) {
       srsi[j] <- 0  # All rows are equal
    } else {
      rankchanges <- sum(diff(ranks) != 0)
       srsi[j] <- rankchanges / (n - 1)
    }
  }
  
  stability_table <- as.data.frame(rbind(SD=sdev, CRV=crv, SD2=sdev2, SRSI= srsi, RSI=stabscores))
  stability_table <- round(stability_table, 2)
  colnames(stability_table) <- colnames(rankmat)
  
  # Weight Sensitivity Analysis (According to reference method)
  sensscores <- numeric(m)
  names(sensscores) <- colnames(rankmat)
  
  for (j in 1:m) {
    rankshifts <- numeric(n)  
    
    for (i in 1:n) {
      rankdiff <- abs(rankmat[i, j] - rankmat[1, j])
      rankshifts[i] <- rankdiff
    }
    sensscores[j] <- mean(rankshifts)
  }
  sensscores <- round(sensscores, 2)
  
  sensscores2 <- numeric(ncol(rankmat))
  for (j in seq_len(ncol(rankmat))) {
    d <- dist(rankmat[, j], method = "manhattan")   
    sensscores2[j] <- round(mean(d), 2)
  }
  names(sensscores2) <- colnames(rankmat)

  
  results <- list(stabtable=stability_table, sensscores=sensscores, sensscores2=sensscores2)
  return(results)
}