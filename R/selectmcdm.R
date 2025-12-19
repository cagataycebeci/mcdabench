selectmcdm <- function (rankmat, dmatrix, bcvec, weights, mval = NULL, 
  mvmet = "percent", mvpars = list(perc = 0.25, rule = c(2, 4, 5), k = 5, r = 3), 
  verbose = TRUE) 
{
  if (!is.matrix(rankmat) && !is.data.frame(rankmat)) stop("Error: 'rankmat' must be a matrix or data.frame.")
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) stop("Error: 'dmatrix' must be a matrix or data.frame.")
  if (nrow(dmatrix) != ncol(rankmat)) stop("Error: The number of rows in 'dmatrix' must match the number of columns in 'rankmat'.")
  if (ncol(dmatrix) != length(bcvec) || ncol(dmatrix) != length(weights)) stop("Error: Lengths of 'bcvec' and 'weights' must match number of criteria.")
  if (!all(bcvec %in% c(-1, 1))) stop("Error: 'bcvec' must contain only -1 or 1.")

  n <- ncol(rankmat)
  
  if (is.null(mval)) {
    if (mvmet == "percent") {
      percvalue <- mvpars$perc
      if (is.null(percvalue) || !is.numeric(percvalue) || percvalue <= 0 || percvalue > 1) 
        stop("Error: 'mvpars$perc' must be between 0 and 1.")
      mval <- max(2, round(n * percvalue))
    } else if (mvmet == "entropy") {
      freq <- table(as.vector(rankmat))
      p <- as.numeric(freq) / length(rankmat)
      entropy_actual <- -sum(p * log(p, 2))
      max_entropy <- ifelse(n <= 1, 0, log(n, 2))
      entropy_norm <- ifelse(max_entropy == 0, 0, entropy_actual / max_entropy)
      mval <- min(n, max(2, round(1 + (n - 1) * entropy_norm)))
    } else if (mvmet == "rule") {
      rule_values <- mvpars$rule
      if (length(rule_values) != 3) stop("Error: 'mvpars$rule' must have length 3.")
      mval <- if (n < 5) rule_values[1] else if (n <= 10) rule_values[2] else rule_values[3]
    } else if (mvmet == "topk") {
      if (mvpars$k > ncol(rankmat)) stop("Parameter k exceeds number of alternatives.")
      rankcounts <- colSums(rankmat <= mvpars$k)
      if (verbose) {
        message("Frequency of alternatives in top k positions:")
        print(rankcounts)
      }
      mval <- max(2, sum(rankcounts >= mvpars$r))
    } else stop("Unsupported 'mvmet'.")
    
    if (verbose) message("Dynamic mval calculated using method '", mvmet, "': ", mval)
  }
  
  if (mval <= 0 || mval > n) stop("Invalid 'mval'.")

  if (verbose) {
    message("* Performance-Based Comparison (PBC) started ...")
    message("  Determined number of top alternatives (m): ", mval)
    message("  Total number of alternatives (n): ", n)
    message("  Number of MCDA techniques used: ", nrow(rankmat), "\n")
  }

  convcritdata <- as.matrix(dmatrix)
  if (verbose) message("Transforming criteria based on type.")
  
  for (j in 1:ncol(convcritdata)) {
    if (bcvec[j] == -1) {
      min_val <- min(convcritdata[, j])
      convcritdata[, j] <- convcritdata[, j] + min_val
      if (verbose) message("  Criterion '", colnames(dmatrix)[j], "' (Cost) transformed. Added ", min_val)
    } else {
      if (verbose) message("  Criterion '", colnames(dmatrix)[j], "' (Benefit) unchanged.")
    }
  }

  weimax <- max(weights)
  highweicritidx <- which(weights == weimax)
  
  if (length(highweicritidx) > 1) {
    if (verbose) message("Step 4: Multiple highest-weighted criteria found. Composite criterion created.")
    compcrit <- rowSums(convcritdata[, highweicritidx, drop = FALSE])
  } else {
    if (verbose) message("Step 4: Single highest-weighted criterion found.")
    compcrit <- convcritdata[, highweicritidx]
  }

  names(compcrit) <- rownames(dmatrix)

  if (verbose) {
    message("  Composite Criterion Values:")
    print(compcrit)
  }

  mcdmcumvals <- list()
  if (verbose) message("Step 5: Calculating cumulative values.")
  
  for (i in 1:nrow(rankmat)) {
    mcda_name <- rownames(rankmat)[i]
    orderedidx <- order(as.numeric(rankmat[i, ]))
    sortcompvals <- compcrit[orderedidx]
    mcdmcumvals[[mcda_name]] <- cumsum(sortcompvals)
  }

  mcumulval <- setNames(sapply(mcdmcumvals, function(x) x[mval]), names(mcdmcumvals))
  
  suggestedmcdm <- names(which(mcumulval == max(mcumulval, na.rm = TRUE)))

  return(list(m = mval,
              suggestedmcdm = suggestedmcdm,
              highcumval = max(mcumulval, na.rm = TRUE),
              all_cumulvalm = mcumulval))
}
