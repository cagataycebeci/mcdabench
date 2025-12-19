\name{rov}
\alias{rov}
\encoding{UTF-8}
\title{ROV: Range of Value}
\description{
  Applies the Rank of Value (ROV) method on a decision matrix using specified weights and a benefit-cost vector.
}
\usage{
rov(dmatrix, bcvec, weights, normethod = "maxmin", tiesmethod="average")
}
\arguments{
  \item{dmatrix}{The decision matrix. Rows represent alternatives, and columns represent criteria.}
  \item{bcvec}{A vector indicating the type of each criterion. Use 1 for benefit criteria and -1 for cost criteria.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"} to apply equal weighting for criteria. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"maxmin"}. See \code{\link{calcnormal}} for details. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\value{
  A list includes the following results:
  \item{u_i}{The utility score obtained for each alternative.}
  \item{rank}{The ranking of the alternatives based on their utility scores. Alternatives with the same score receive the average rank.}
}
\details{
  The \code{rov} function implements the Rank Order Centroid (ROV) method, a multi-criteria decision-making (MCDM) technique used to rank alternatives. The function takes a decision matrix, a benefit-cost vector, and criterion weights as input. The decision matrix is first normalized using the specified normalization method. Then, the criterion weights are either provided directly or calculated using a specified method. Finally, a utility score is calculated for each alternative, and the alternatives are ranked based on these scores.
}
\references{
Madic, M., Radovanovic, M., & Manic, M. (2016). Application of the ROV method for the selection of cutting fluids. Decision Science Letters, 5(2), 245-254.
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\author{
  Cagatay Cebeci
}
\examples{
# Sample decision matrix
dmat <- matrix(c(
  84, 10, 5, 20,  
  66, 18, 8, 15,
  74, 12, 6, 25,
  90, 22, 9, 18, 
  68, 18, 7, 15  
), nrow = 5, byrow = TRUE)
rownames(dmat) <- c("A1", "A2", "A3", "A4", "A5")
colnames(dmat) <- c("C1", "C2", "C3", "C4")

# Benefit-cost vector (1: benefit, -1: cost)
bc <- c(1, -1, 1, -1)

# User-defined weights
uw <- c(0.1, 0.4, 0.3, 0.2)
  
# Apply ROV using user-defined weights
resrov_1 <- rov(dmatrix=dmat, bcvec=bc, weights=uw)
print(resrov_1)

# Apply ROV using equal weights
equw <- calcweights(dmatrix=dmat, bcvec=bc, type="equal", normethod="maxmin")
resrov_2 <- rov(dmatrix=dmat, bcvec=bc, weights = equw)
print(resrov_2)
}

\keyword{ multivariate }
\keyword{ decision-making }