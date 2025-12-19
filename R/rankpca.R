rankpca <- function(rankmat, biplot = TRUE, reverse=FALSE) {
  stopifnot(is.matrix(rankmat) || is.data.frame(rankmat))
  if (anyNA(rankmat)) stop("rankmat cannot contain NA")
  if (is.null(rownames(rankmat)) || ncol(rankmat) < 1)
    stop("rankmat needs row names and at least one column")

  if(reverse)
    mat <- (ncol(rankmat) + 1) - as.matrix(rankmat)
  else
    mat <- rankmat

  keep <- apply(mat, 2, sd) > 0
  if (!all(keep)) {
    warning("Constant columns removed: ",
            paste(colnames(mat)[!keep], collapse = ", "))
    mat <- mat[, keep, drop = FALSE]
  }
  if (ncol(mat) < 2) {
    warning("After removing constant columns fewer than 2 variables remain; PCA skipped.")
    return(NULL)
  }

  pca <- prcomp(mat, scale. = TRUE, center = TRUE)

  results <- list(
    respcas   = pca,
    ascores   = as.data.frame(pca$x),
    critloads = as.data.frame(pca$rotation),
    expvar    = summary(pca)$importance[2, ],
    cumvar    = summary(pca)$importance[3, ],
    intprtpc1 = if (sign(mean(pca$rotation[, 1])) <= 0)
                  "Higher PC1 scores typically indicate alternatives with generally higher ranks (worse performance) across criteria. Alternatives on the left of the biplot are generally better, those on the right are generally worse."
                else
                  "Higher PC1 scores typically indicate alternatives with generally lower ranks (better performance) across criteria. Alternatives on the left of the biplot are generally worse, those on the right are generally better.",
    intprtpc2 = "PC2 scores close to zero (near the horizontal axis) indicate alternatives with uniform ranking behavior across criteria. Scores significantly greater or lower than zero (further from the horizontal axis) suggest non-uniform (variable) performance across different criteria."
  )

  if (biplot){
    p <- factoextra::fviz_pca_biplot(
      pca, repel = TRUE,
      col.ind = "darkblue", col.var = "darkred",
      geom = c("point", "text"), labelsize = 4,
      title = "PCA biplot")
    print(p)
  }
  return(results)
}
