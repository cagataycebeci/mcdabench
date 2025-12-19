\name{aras}
\alias{aras}
\title{ARAS: Additive Ratio Assessment Method}
\description{
The ARAS method, developed by Zavadskas and Turskis in 2010, aims to evaluate and rank alternatives based on given criteria using normalization, weighting, and utility degree computation. Due to its simple calculation procedures and clear logic, ARAS is widely used in multi-criteria decision-making (MCDA) problems.
}
\usage{
aras(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix, where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion as a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The ARAS method extends the decision matrix with ESP values, calcnormals the scores, applies weights, calculates utility values (\code{S}) 
and relative utilities (\code{K}) for ranking alternatives.
}
\value{
A list containing the following components:
  \item{Si}{The total utility values of the alternatives.}
  \item{Ki}{The relative utility values for ranking the alternatives.}
  \item{rank}{Ranking the alternatives according to the relative utility values.}
}
\references{
Zavadskas, E. K., & Turskis, Z. (2010). A new additive ratio assessment (ARAS) method in multicriteria decision-making. \emph{Technological and Economic Development of Economy},
\emph{16}(2), 159-172. <doi:10.3846/tede.2010.10>.
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

# Apply ARAS with user-defined weights
resaras_1 <- aras(dmatrix=dmat, bcvec=bc, weights=uw)
print(resaras_1)

# Apply ARAS with CRITIC weights
criwei <- calcweights(dmatrix=dmat, bcvec=bc, type="critic", normethod="ratio")

resaras_2 <- aras(dmatrix=dmat, bcvec=bc, weights=criwei)
print(resaras_2$Ki)    # Relative utilities
print(resaras_2$rank)  # Ranks of relative utilities
}
