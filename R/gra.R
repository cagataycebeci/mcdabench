gra <- function(dmatrix, bcvec, weights, idesol = NULL, grdmethod="sum", rho=0.5, tiesmethod="average") {
  # Input checks
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("The sum of the weights must be approximately 1.")
  }
  if (rho < 0 || rho > 1) {
    stop("The rho parameter must be between 0 and 1. Try a value like 0.5 for a typical setting.")
  }
  if (!grdmethod %in% c("sum", "mean")) {
    stop('Invalid method for calculating Grey Relation Degrees. Choose either "sum" or "mean".')
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)
  
  nmatrix <- calcnormal(dmatrix, bcvec, type="maxmin")

  # Determining the Ideal Solution
  if (is.null(idesol)) {
    idesol <- numeric(m)
    for (j in 1:m) {
      if (bcvec[j] == 1) {
        idesol[j] <- max(nmatrix[, j])
      } else {
        idesol[j] <- min(nmatrix[, j])
      }
    }
  } else {
    if (length(idesol) != m) {
      stop("The length of ideal solutions (idesol) must be equal to the number of columns in dmatrix.")
    }
  }

  # Calculating Grey Relational Coefficients
  greycoeffs <- matrix(0, nrow = n, ncol = m)
  min_difference <- apply(abs(idesol - nmatrix), 2, min) + 1e-6
  max_difference <- apply(abs(idesol - nmatrix), 2, max) + 1e-6

  for (i in 1:n) {
    for (j in 1:m) {
      difference <- abs(idesol[j] - nmatrix[i, j])
      greycoeffs[i, j] <- (min_difference[j] + rho * max_difference[j]) / 
         (difference + rho * max_difference[j])
    }
  }

  # Calculating Grey Relation Degrees based on method
  if (grdmethod == "sum") {
    greyreldegrees <- rowSums(greycoeffs * weights)
  } else {
    greyreldegrees <- rowMeans(greycoeffs * weights)
  }

  ranking <- rank(-greyreldegrees, ties.method = tiesmethod)
  names(ranking) <- rownames(dmatrix)
  
  ranked_alternatives <- data.frame(Alternative = rownames(dmatrix),
          Grey_Relation_Degree = greyreldegrees,
          Rank = ranking)
  ranked_alternatives <- ranked_alternatives[order(ranked_alternatives$Rank), ]

  results <- list(
    Ideal_Solution = idesol,
    Grey_Relation_Coef = greycoeffs,
    Grey_Relation_Deg = greyreldegrees,
    Ranking = ranked_alternatives,
    rank = ranking
  )
  
  return(results)
}
