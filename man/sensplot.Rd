\name{sensplot}
\alias{sensplot}
\encoding{UTF-8}
\title{Plot Weight Sensitivity Analysis Results}
\description{
Creates a bar or pie chart to visualize the frequency of ranking patterns resulting from a sensitivity analysis performed by \code{\link{weisana}}.
}
\usage{
sensplot(senstable, topn = 10, colpal = NULL, type = "bar", mtitle = NULL)
}
\arguments{
  \item{senstable}{A data frame containing ranking patterns and their frequencies, typically the \code{sensitivity_table} output from \code{weisana()}. Must include \code{Pattern} and \code{Percent} columns.}
  \item{topn}{Number of top-ranked patterns to include in the plot. Defaults to \code{10}.}
  \item{colpal}{Optional color palette for the bars or pie slices. If \code{NULL}, a \code{rainbow()} palette is used.}
  \item{type}{Type of plot. Choose between \code{"bar"} (default) or \code{"pie"}.}
  \item{mtitle}{Main title of the plot. If \code{NULL}, defaults to \code{"Sensitivity Analysis"}.}
}
\details{
This function is a companion to \code{weisana()} and is designed to offer a quick graphical summary of how frequently different ranking outcomes appear across weight scenarios. The bar plot also includes percentage labels above the bars for easy interpretation.
}
\examples{
# Simulated example
set.seed(123)
results <- data.frame(
  Pattern = c("2,1,2,4", "2,1,3,4", "1,3,2,4", "1,3,4,2"),
  Count = c(50, 30, 15, 5),
  Percent = c(50, 15, 20, 5),
  stringsAsFactors = FALSE
)

sensplot(results, topn = 5, type = "bar")
sensplot(results, type = "pie", mtitle = "Ranking Pattern Distribution")
}
\seealso{
\code{\link{weisana}} for generating the sensitivity table input.
}
\author{Cagatay Cebeci}
