\name{mabac}
\alias{mabac}
\encoding{UTF-8}
\title{MABAC: Multi-Attributive Border Approximation Area Comparison}
\description{
This function is an implementation of the Multi-Attributive Border Approximation Area Comparison (MABAC) method to evaluate and rank alternatives based on given criteria through normalization, weighting, and border approximation area calculations.
}
\usage{
mabac(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"}. See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The MABAC method involves normalizing the decision matrix, applying weights to criteria, and calculating the border approximation area. The border approximation area represents the average value of the alternatives for each criterion. A distance matrix is computed by subtracting the border approximation area from the weighted normalized matrix, and the alternatives are ranked based on their scores.
}
\value{
A list containing the following components:
  \item{G}{The border approximation area for each criterion.}
  \item{Q}{The distance matrix for each alternative and criterion.}
  \item{score}{The final scores for each alternative.}
  \item{rank}{The ranks of final scores of alternatives.}
}
\references{
Pamučar, D., & Ćirović, G. (2015). The selection of transport and handling resources in logistics centers using Multi-Attributive Border Approximation area Comparison (MABAC). 
\emph{Expert Systems with Applications, 42}(6), 3016-3028. <doi:10.1016/j.eswa.2014.11.057>
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

# Apply MABAC with user-defined weights
resmabac_1 <- mabac(dmatrix=dmat, bcvec=bc, weights=uw)
print(resmabac_1) 

# Apply MABAC with MEREC weighting method
mw <- calcweights(dmat, bcvec=bc, type="merec", normethod="maxmin")
resmabac_2 <- mabac(dmatrix=dmat, bcvec=bc, weights=mw)
print(resmabac_2$score) # Relative utilities
print(resmabac_2$rank)  # Ranks of relative utilities
}
