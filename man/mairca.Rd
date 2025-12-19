\name{mairca}
\alias{mairca}
\encoding{UTF-8}
\title{MAIRCA: Multi-Attribute Ideal-Real Comparative Analysis}
\description{
  The MAIRCA is a multi-criteria decision-making method that normalizes the decision matrix and ranks alternatives based on specified criterion weights.
}
\usage{
  mairca(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"}. See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\value{
  Returns a list containing:
  \item{Q}{Final values representing the ranking of alternatives.}
  \item{rank}{Ranking of the alternatives.}
}
\references{
Gigović, L., Pamučar, D., Bajić, Z., & Milićević, M. (2016). The combination of expert judgment and GIS-MAIRCA analysis for the selection of sites for ammunition depots.
 \emph{Sustainability, 8}(4), 372. <doi:10.1080/1331677X.2018.1506706>
}
\seealso{
  \code{\link{calcweights}}
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
  
# Apply MAIRCA with user-defined weights
resmairca_1 <- mairca(dmatrix=dmat, bcvec=bc, weights = uw)
print(resmairca_1)

# Apply MAIRCA with Critic weights
cw <- calcweights(dmat, bcvec=bc, type="critic")
resmairca_2 <- mairca(dmatrix=dmat, bcvec=bc, weights = cw)
print(resmairca_2)
}