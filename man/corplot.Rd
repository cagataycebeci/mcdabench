\name{corplot}
\alias{corplot}
\encoding{UTF-8}
\title{Correlation Matrix Visualization}
\description{
Visualizes a correlation or similarity matrix using heatmap-style plots, allowing customization of colors, shapes, axis labels, and text annotations.
}
\usage{
corplot(cormat, labsize = 10, fsize = 3, fcolor = "black", 
   colpal = c("red", "white", "dodgerblue"), title = "Correlation Matrix", 
   xlab = "Coef", ylab = "Coef", shape = "square", displabels = TRUE)
}
\arguments{
  \item{cormat}{A matrix representing correlations or similarities.}
  \item{labsize}{Numeric, base size of plot text and labels. Default is 10.}
  \item{fsize}{Numeric, font size for cell text annotations. Default is 3.}
  \item{fcolor}{Character, color of text annotations. Default is "black".}
  \item{colpal}{Character vector, specifying three colors for low, mid, and high values. Default is c("red", "white", "blue").}
  \item{title}{Character, title of the plot. Default is "Correlation Matrix".}
  \item{xlab}{Character, label for the X-axis. Default is "Coef".}
  \item{ylab}{Character, label for the Y-axis. Default is "Coef".}
  \item{shape}{Character, specifying either "circle" or "square" for visualization. Default is "square".}
  \item{displabels}{Logical, whether to display correlation values inside cells. Default is TRUE.}
}
\examples{
# Sample correlation matrix
set.seed(123)
df <- data.frame(
  X1 = rnorm(20), X2 = rnorm(20), X3 = rnorm(20), X4 = rnorm(20),
  X5 = rnorm(20), X6 = rnorm(20), X7 = rnorm(20), X8 = rnorm(20),
  X9 = rnorm(20), X10 = rnorm(20)
)
cormat <- cor(df, method = "spearman")

corplot(cormat, labsize = 12, fsize = 4, fcolor = "black", 
        title = "Spearman Correlation Matrix", xlab = "Methods", ylab = "Methods",
        shape = "square", displabels = TRUE)
}
