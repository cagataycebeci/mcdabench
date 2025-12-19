\name{copras}
\alias{copras}
\encoding{UTF-8}
\title{COPRAS: Complex Proportional Assessment}
\description{
 Applies the COPRAS method for decision makings.
}
\usage{
copras(dmatrix, bcvec, weights, normethod = "sum", tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{A vector indicating the criteria type: 1 for benefit criteria, -1 for cost criteria.}
  \item{weights}{A numeric vector representing the relative importance of each criterion, summing to 1.}
  \item{normethod}{Normalization method applied to the decision matrix. The \code{"sum"} normalization
  is default to COPRAS.}
  \item{tiesmethod}{Method to resolve ties in ranking, default is \code{"average"}.}
}
\value{
A list with the following components:
  \item{weights}{The weights applied to each criterion.}
  \item{splus}{The positive significance sum for each alternative.}
  \item{sminus}{The negative significance sum for each alternative.}
  \item{qi}{The relative significance (Q) values.}
  \item{ui}{Normalized utility values of each alternative.}
  \item{rank}{Ranking of alternatives based on utility values.}
}
\references{
Razavi Hajiagha, S. H., Hashemi, S. S., & Zavadskas, E. K. (2013). A complex proportional assessment method for group decision making in an interval-valued intuitionistic fuzzy environment. 
\emph{Technological and Economic Development of Economy, 19}(1), 22-37. <doi:10.3846/20294913.2012.762953>
}
\examples{
dmatrix <- matrix(c(
  10, 50, 30, 80,
  20, 60, 20, 70,
  30, 40, 40, 60,
  40, 30, 50, 50
), nrow=4, byrow=TRUE)

bcvec <- c(1, 1, -1, 1)
weights <- c(0.3, 0.4, 0.1, 0.2)

rescopras <- copras(dmatrix, bcvec, weights)
print(rescopras)
}

\author{
Cagatay Cebeci
}