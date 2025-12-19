# COmplex PRoportional ASsessment (COPRAS) 
copras <- function(dmatrix, bcvec, weights, normethod="sum", tiesmethod="average") {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(-1, 1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain -1 (cost) or 1 (benefit) values and have the same length as the number of criteria.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("The weights must be a numeric vector with the same length as the number of criteria.")
  }

  # Normalization (Sum Normalization)
  #nmatrix <- sweep(dmatrix, 2, colSums(dmatrix), "/")
  nmatrix <- calcnormal(dmatrix,bcvec, type=normethod)

  # Weighted Normalized Matrix 
  wmatrix <- sweep(nmatrix, 2, weights, "*")

  # Identify indices for benefit and cost criteria
  idxmax <- which(bcvec == 1)  # Benefit criteria indices
  idxmin <- which(bcvec == -1) # Cost criteria indices

  # Calculate S+ and S- Values 
  splus <- rowSums(wmatrix[, idxmax, drop = FALSE])
  sminus <- rowSums(wmatrix[, idxmin, drop = FALSE])

  # Division by zero prevention by replacing zeros with a small epsilon
  sminus[sminus == 0] <- 1e-6

  # Calculate Relative Significance (Q) Values
  if (length(idxmin) > 0) { # If there are cost criteria
    sum_sminus_overall <- sum(sminus)
    sum_inv_sminus <- sum(1 / sminus)
    qi <- splus + (sum_sminus_overall / (sminus * sum_inv_sminus))
  } else {
    qi <- splus
  }

  # Utility Values
  ui <- qi / max(qi)
  ranking <- rank(-ui, ties.method = tiesmethod)

   results <- list(
    weights = weights,
    splus = splus,
    sminus = sminus,
    qi = qi,
    ui = ui,
    rank = ranking
  )

  return(results)
}

