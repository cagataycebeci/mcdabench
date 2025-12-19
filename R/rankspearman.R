rankspearman <- function(rankmat) {
   n <- nrow(rankmat)
  
  # Compute the Spearman correlation matrix for rows
  src <- round(suppressWarnings(cor(t(rankmat), method = "spearman")), 2)  # Transpose for row-wise correlation
  
  # Initialize a matrix to store p-values
  pvals <- matrix(NA, nrow = n, ncol = n)
  colnames(pvals) <- rownames(rankmat)
  rownames(pvals) <- rownames(rankmat)

  # Loop through all row pairs to calculate p-values
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        test <-  suppressWarnings(cor.test(rankmat[i, ], rankmat[j, ], method = "spearman", exact = FALSE))
        pvals[i, j] <- test$p.value
      } else {
        pvals[i, j] <- 1  # Set diagonal values to 1
      }
    }
  }

  # Replace NA values in correlation matrix with 1
  src[is.na(src)] <- 1  
  pvals[is.na(pvals)] <- 1  

  # Create a display matrix for visualization
  displaymat <- matrix("", nrow = n, ncol = n)
  colnames(displaymat) <- rownames(rankmat)
  rownames(displaymat) <- rownames(rankmat)

  # Populate the upper triangle with significance stars, and lower triangle with correlation values
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        displaymat[i, j] <- "1"  # Diagonal elements set to 1
      } else if (i < j) {  # Upper triangle (significance levels)
        if (!is.na(pvals[i, j])) {
          if (pvals[i, j] < 0.001) {
            displaymat[i, j] <- "***"
          } else if (pvals[i, j] < 0.01) {
            displaymat[i, j] <- "**"
          } else if (pvals[i, j] < 0.05) {
            displaymat[i, j] <- "*"
          } else {
            displaymat[i, j] <- "" 
          }
        }
      } else {  # Lower triangle (correlation values)
        displaymat[i, j] <- format(src[i, j], nsmall = 2)  
      }
    }
  }
  displaymat <- as.data.frame(displaymat)
  
  results <- list(
    cormat = src,  # Spearman correlation matrix
    p_values = pvals,  # Matrix of p-values
    displaymat = displaymat  # Formatted visualization matrix
  )

  return(results)
}
