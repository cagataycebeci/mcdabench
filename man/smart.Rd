\name{smart}
\alias{smart}
\encoding{UTF-8}
\title{SMART: Simple Multi-Attribute Rating Technique}
\description{
  This function applies the Simple Multi-Attribute Rating Technique (SMART) for multi-criteria decision analysis. The SMART evaluates alternatives based on weighted aggregation of their values determined by value functions.
}
\usage{
smart(dmatrix, bcvec=NULL, weights, normethod=NULL, valfuncs = NULL, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}. Default is \code{NULL} because SMART uses a different techique.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"} to apply equal weighting for criteria. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{NULL}. See \code{\link{calcnormal}} for details.}
  \item{valfuncs}{An optional list of value functions for each criterion. Each element of the list should be a function that takes a numeric value (criterion performance) and returns a scaled value (typically between 0 and 1). If \code{NULL} (default), linear value functions (MaxMin scaling) are applied to each criterion.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
  The SMART method involves defining value functions for each criterion to scale the performance values and assigning weights to the criteria to reflect their importance. The overall score for each alternative is then calculated as the weighted sum of its values across all criteria. Alternatives are ranked based on their overall scores in descending order.

  If \code{weights} is a numeric vector, these weights are used directly after normalization.

  If \code{valfuncs} is not provided, linear value functions (min-max scaling) are used to normalize the performance values of each criterion to a 0-1 range.
}
\value{
  A list containing the following components:
  \item{Value_Matrix}{A numeric matrix containing the scaled values of each alternative for each criterion based on the value functions.}
  \item{Weighted_Value_Matrix}{A numeric matrix containing the weighted scaled values of each alternative for each criterion.}
  \item{Overall_Scores}{A numeric vector of the overall scores for each alternative.}
  \item{Ranking}{A data frame with the ranked alternatives, their SMART scores, and their ranks.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower rank is better).}
}
\references{
Olson, D.L. (1996). Smart. In: \emph{Decision Aids for Selection Problems}. Springer Series in Operations Research. Springer, New York, NY. <doi:10.1007/978-1-4612-3982-6_4>
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
ressmart_1 <- smart(dmatrix=dmat, bcvec=bc, weights=uw)
print(ressmart_1)

# Criteria weights with std. dev. method
sdw <- calcweights(dmat, bcvec=bc, type="sdev")
ressmart_2 <- smart(dmatrix=dmat, bcvec=bc, weights = sdw)
print(ressmart_2)

# Custom value functions
vf <- list(
  function(x) 1 - (x - min(dmat[, 1])) / (max(dmat[, 1]) - min(dmat[, 1])),
  function(x) (x - min(dmat[, 2])) / (max(dmat[, 2]) - min(dmat[, 2])),
  function(x) 1 - (x - min(dmat[, 1])) / (max(dmat[, 1]) - min(dmat[, 1])),
  function(x) (x - min(dmat[, 2])) / (max(dmat[, 2]) - min(dmat[, 2]))
)
ressmart_3 <- smart(dmatrix=dmat, bcvec=bc, weights=uw, valfuncs = vf)
print(ressmart_3)
}