\name{wsm}
\alias{wsm}
\encoding{UTF-8}
\title{WSM: Weighted Sum Method / SAW: Simple Additive Weighting}
\description{
Implementation of the Weighted Sum Method (WSM) or Simple Addiitive Weighting (SAW) to evaluate and rank alternatives by aggregating weighted normalized scores for each criterion.
}
\usage{
wsm(dmatrix, bcvec, weights, normethod = NULL, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{
    Decision matrix where rows represent alternatives and columns represent criteria.
  }
  \item{bcvec}{
    Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).
  }
  \item{weights}{A numeric vector of user-defined weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{NULL}. See \code{\link{calcnormal}} for details. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
Simple Additive Weighting method emerged from the early developments in multi-attribute utility theory and decision analysis in the 1960s, with significant contributions from researchers like Ward Edwards and David E. Bell. The method was then further popularized and applied in various contexts within the MCDM field. The WSM evaluates alternatives by normalizing the decision matrix, applying weights to criteria using multiplication, and aggregating the weighted scores (\code{p}) for each alternative. This straightforward method provides a clear ranking of alternatives based on their overall performance.
}
\value{
A list containing the following components:
  \item{nmatrix}{The normalized decision matrix.}
  \item{weighted_matrix}{The decision matrix after applying weights.}
  \item{p}{The aggregated weighted scores for each alternative.}
  \item{rank}{The ranks of aggregated weighted scores for alternatives.}
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\references{
Zavadskas, E. K., Vilutiene, T., Turskis, Z., & Šaparauskas, J. (2014). Multi-criteria analysis of Projects' performance in construction. 
 \emph{Archives of Civil and Mechanical Engineering, 14}, 114-121. <doi:10.1016/j.acme.2013.07.006>

Panjaitan, M. I. (2019). Simple Additive Weighting (SAW) method in determining beneficiaries of foundation benefits. \emph{Login, 13}(1), 19-25. <doi:10.24224/login.v13i1.22>
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

reswsm_1 <- wsm(dmatrix=dmat, bcvec=bc, weights=uw)
print(reswsm_1)

reswsm_2 <- wsm(dmatrix=dmat, bcvec=bc, weights=uw, normethod = "sum")
print(reswsm_2)

reswsm_3 <- wsm(dmatrix=dmat, bcvec=bc, weights=uw, normethod = "maxmin")
print(reswsm_3)
}
