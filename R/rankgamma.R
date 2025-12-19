# Calculates Goodman & Kruskal's Gamma statistic between two rankings
rankgamma <- function(rankmat) {
  n <- nrow(rankmat)
  gmatrix <- matrix(NA, nrow = n, ncol = n)
  rownames(gmatrix) <- colnames(gmatrix) <- rownames(rankmat)

  gkg <- function(x, y) {
    concordant <- 0
    discordant <- 0
    m <- length(x)
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        dx <- sign(x[i] - x[j])
        dy <- sign(y[i] - y[j])
        if (dx != 0 && dy != 0) {
          if (dx == dy) concordant <- concordant + 1
          else discordant <- discordant + 1
        }
      }
    }
    denom <- concordant + discordant
    if (denom == 0) return(NA)
    return((concordant - discordant) / denom)
  }

  for (i in 1:n) {
    for (j in 1:n) {
      gmatrix[i, j] <- gkg(rankmat[i, ], rankmat[j, ])
    }
  }

  return(list(rankmat = rankmat, gammamat = gmatrix))
}




