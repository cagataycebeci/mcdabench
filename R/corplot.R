# Correlation plot
utils::globalVariables(c("value", "Var1", "Var2"))
corplot <- function(cormat, labsize = 10, fsize = 3, fcolor = "black", 
   colpal = c("red", "white", "dodgerblue"), title = "Correlation Matrix", 
   xlab = "Coef", ylab = "Coef",  shape = "square", displabels = TRUE) {
  # Input checks
  if (!is.matrix(cormat)) {
    stop("Input data (cormat) must be a matrix representing a correlation or similarity matrix.")
  }
  if (!is.character(shape) || !shape %in% c("circle", "square")) {
    stop("shape must be one of 'circle' or 'square'.")
  }
  if (!is.logical(displabels)) {
    stop("displabels must be TRUE or FALSE.")
  }

  long_df <- as.data.frame(expand.grid(Var1 = rownames(cormat), Var2 = colnames(cormat)))
  long_df$value <- as.vector(cormat)

  if (shape == "circle") {
     geom_shape <- ggplot2::geom_point(ggplot2::aes(size = abs(value), fill = value), 
          shape = 21, color = "black", stroke = 0.2)
  } else {
     geom_shape <- ggplot2::geom_tile(color = "black", size = 0.2)
  }

   p <- ggplot2::ggplot(long_df, ggplot2::aes(Var1, Var2, fill = value)) +
    geom_shape +
    ggplot2::scale_fill_gradient2(low = colpal[1], mid = colpal[2], high = colpal[3], midpoint = 0) +
    ggplot2::scale_size_continuous(range = c(2, 10)) + 
    ggplot2::labs(x = xlab, y = ylab, title = title) + 
    ggplot2::theme_minimal(base_size = labsize) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  if (displabels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), color = fcolor, size = fsize)
  }

  print(p)
}
