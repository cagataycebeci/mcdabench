\name{ram}
\alias{ram}
\encoding{UTF-8}
\title{RAM: Root Assessment Method}
\description{
Implements the Root Assessment Method (RAM) to evaluate and rank alternatives based on normalized scores, weighted contributions, and root-based calculations.
}
\usage{
ram(dmatrix, bcvec, weights, normethod = "sum", tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"sum"}. See \code{\link{calcnormal}} for details. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The RAM method evaluates alternatives by normalizing the decision matrix, applying weights to criteria, and calculating scores based on the root formula that integrates positive and negative contributions. The final ranking is derived from the computed scores (\code{ri}), which reflect the balance between benefit and cost criteria.
}
\value{
A list containing the following components:
  \item{nmatrix}{The normalized decision matrix.}
  \item{weighted_matrix}{The decision matrix after applying weights.}
  \item{Spi}{The summed positive contributions for each alternative.}
  \item{Smi}{The summed negative contributions for each alternative.}
  \item{ri}{The final root-based scores for each alternative.}
  \item{rank}{The ranks of final root-based scores for alternatives.}
}
\references{
Sotoudeh-Anvari, A. (2023). Root Assessment Method (RAM): A novel multi-criteria decision making method and its applications in sustainability challenges. \emph{Journal of Cleaner Production, 423}, 138695.
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
  
# Apply RAM with user-defined weights
resram_1 <- ram(dmatrix=dmat, bcvec=bc, weights=uw)
print(resram_1)

# Apply RAM with equal weights
equw <- calcweights(dmat, bcvec=bc, type="equal")
resram_2 <- ram(dmatrix=dmat, bcvec=bc, weights=equw)
print(resram_2)
}
