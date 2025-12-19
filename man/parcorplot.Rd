\name{parcorplot}
\alias{parcorplot}
\title{Parallel Coordinate Plot}
\description{
  This function creates a parallel coordinate plot using the provided data frame or matrix.
  Each row represents a different method, and each line in the plot corresponds to the values
  for that method.
}
\usage{
parcorplot(x, xl = NULL, yl = NULL, lt = NULL, colpal = NULL)
}
\arguments{
  \item{x}{A data frame or matrix where the row names are the names of MCDM methods,
    and the columns contain the corresponding rank values for alternatives computed
    by the methods.}
  \item{xl}{Label for the x-axis. Default is "x".}
  \item{yl}{Label for the y-axis. Default is "y".}
  \item{lt}{Title for the legend. Default is "x".}
  \item{colpal}{A vector of colors for the lines. If NULL (default),
    a predefined color palette or \code{topo.colors} is used.}
}
\value{
  No return value, but a plot displaying the parallel coordinates for each method.
}
\author{
  Cagatay Cebeci
}
\examples{
# Sample rank matrix
rankmat <- data.frame(
  ALT1 = c(3.0, 3.0, 2.5, 3.0, 3.0, 3.0, 2.0, 1.0),
  ALT2 = c(1, 1, 1, 1, 1, 1, 1, 2),
  ALT3 = c(2.0, 2.0, 2.5, 2.0, 2.0, 2.0, 3.0, 3.0)
)
rownames(rankmat) <- c(
  "ARAS", "EDAS", "ELECTRE4", "MABAC", "MARCOS", "MEGAN", "TOPSIS", "VIKOR"
)

parcorplot(rankmat, xl = "Alternatives", yl = "Ranks", lt = "MCDA Methods")
}