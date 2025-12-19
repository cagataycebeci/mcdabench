rankwilcox <- function(rankmat, alpha = 0.05, padjmethod = "none") {
  n <- nrow(rankmat)
  rnames <- rownames(rankmat)
  pvalmat <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))
  teststatmat <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))
  pvals <- numeric()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      group1 <- rankmat[i, ]
      group2 <- rankmat[j, ]
      wilcox_test <- wilcox.test(group1, group2, paired = TRUE, exact = FALSE)
      pval <- wilcox_test$p.value
      test_stat <- round(wilcox_test$statistic, 3)

      pvalmat[i, j] <- pval
      pvalmat[j, i] <- pval
      teststatmat[i, j] <- test_stat
      teststatmat[j, i] <- test_stat

      pvals <- c(pvals, pval)
    }
  }

  adjpvals <- p.adjust(pvals, method = padjmethod)
  adjpvalmat <- matrix(NA, nrow = n, ncol = n, dimnames = list(rnames, rnames))
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      adjusted_pval <- adjpvals[k]
      adjpvalmat[i, j] <- adjusted_pval
      adjpvalmat[j, i] <- adjusted_pval
      k <- k + 1
    }
  }

  displaymat <- matrix("", nrow = n, ncol = n, dimnames = list(rnames, rnames))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        adjusted_pval <- adjpvalmat[i, j]
        displaymat[i, j] <- format(round(adjusted_pval, 2), nsmall = 2)
      } else if (i > j) {
        displaymat[i, j] <- format(round(teststatmat[i, j], 2), nsmall = 2)
      } else {
        displaymat[i, j] <- ""
      }
    }
  }

  displaymat <- as.data.frame(displaymat)

  return(list(
    alpha = alpha,
    p_adjmethod = padjmethod,
    pvalmat = pvalmat,
    adjpvalmat = adjpvalmat,
    displaymat = displaymat
  ))
}