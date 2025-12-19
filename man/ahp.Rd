\name{ahp}
\alias{ahp}
\encoding{UTF-8}
\title{AHP: Analytic Hierarchy Process}
\description{
  Computes criteria weights using the Analytic Hierarchy Process (AHP). Supports multiple weighting methods and consistency ratio calculation.
}
\usage{
ahp(pmatrix, method = "maxeig", weights = NULL, tiesmethod = "average")
}
\arguments{
  \item{pmatrix}{A square reciprocal comparison matrix.}
  \item{method}{Weight calculation method. Options: \code{"maxeig"}, \code{"mean"}, \code{"geometric"}. Default is \code{"maxeig"}.}
  \item{weights}{Optional predefined weights. If provided, AHP calculations will be skipped.}
  \item{tiesmethod}{Ranking method for handling ties in weight values. Default is \code{"average"}.}
}
\details{
  The AHP method calculates weights based on a comparison matrix.
  If \code{weights} are provided, direct ranking is performed without AHP computation.
  Otherwise, weights are computed using \code{method}, with \code{"maxeig"} using the principal eigenvector.
  The consistency ratio (CR) is also calculated to assess decision matrix reliability.
}
\value{
  A list with the following components:
  \item{weights}{Computed weights for criteria.}
  \item{rc}{Consistency ratio (CR) for matrix reliability.}
  \item{rank}{Ranking of criteria based on weights.}
}
\references{
Saaty, T. L. (1990). How to make a decision: the analytic hierarchy process. \emph{European Journal of Operational Research, 48}(1), 9-26. <doi:10.1016/0377-2217(90)90057-I>
}
\examples{
dataset <- matrix(c(
  1, 6, 2, 5, 7,   
  1/6, 1, 4, 8, 3,  
  1/2, 1/4, 1, 9, 6,  
  1/5, 1/8, 1/9, 1, 2,  
  1/7, 1/3, 1/6, 1/2, 1  
), nrow=5, byrow=TRUE)

resahp_1 <- ahp(dataset, method="maxeig")
print(resahp_1)

resahp_2 <- ahp(dataset, method="mean")
print(resahp_2)
}
\author{Cagatay Cebeci}
