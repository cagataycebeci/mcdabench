\name{boxplotmcda}
\alias{boxplotmcda}
\title{MCDA Boxplot}
\description{
  A function to create boxplots for a Multi-criteria Decision Analysis (MCDA) data matrices with customizable labels, colors, and additional arguments.
}
\usage{
boxplotmcda(datamatrix, labels = NULL, colors = NULL, ax = NULL, mt = NULL)
}
\arguments{
  \item{datamatrix}{A numeric data matrix to generate boxplots.}
  \item{labels}{Optional character vector for labeling the x-axis. Default is \code{NULL}, which assigns sequential labels such as \code{"A1"}, \code{"A2"}, etc.}
  \item{colors}{Optional vector specifying colors for each box. Default is \code{NULL}, which uses a rainbow color scheme.}
  \item{ax}{Optional graphics device. Default is \code{NULL}, which uses the current device.}
  \item{mt}{Optional main title for the plot.}
}
\details{
  This function generates boxplots for the columns of the given data matrix. It allows customization through labels, colors, and additional boxplot arguments. 
  The labels are displayed on the x-axis, and a grid is added for better visual clarity.
}
\author{
  Cagatay Cebeci
}
\examples{
set.seed(123)
datamatrix <- matrix(runif(50, 0, 1), nrow = 10, ncol = 5)
bcvec <- c(-1, 1, 1, -1)

# Data boxplot with default labels and colors
boxplotmcda(datamatrix)

# Data boxplot with custom labels and colors
labels <- c("Criteria 1", "Criteria 2", "Criteria 3", "Criteria 4", "Criteria 5")
colors <- c("red", "blue", "green", "purple", "orange")
boxplotmcda(datamatrix, labels = labels, colors = colors)

# Normalized data boxplot with with default labels and colors
nmatrix <- calcnormal(X=datamatrix, bcvec=bcvec, type="maxmin")
boxplotmcda(nmatrix, mt="MaxMin Normalized Matrix")
}