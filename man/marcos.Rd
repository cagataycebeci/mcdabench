\name{marcos}
\alias{marcos}
\encoding{UTF-8}
\title{MARCOS: Measurement of Alternatives and Ranking According to Compromise Solution}
\description{
Implements the MARCOS method for MCDM by evaluating alternatives based on their relationships to ideal and anti-ideal solutions.
}
\usage{
marcos(dmatrix, bcvec, weights, normethod = "ratio", tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"}. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"ratio"} for MARCOS. See \code{\link{calcnormal}} for details. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The MARCOS method (Stević et al, 2020) constructs an extended decision matrix by adding ideal and anti-ideal solutions. The matrix is normalized, and weights are applied to each criterion. 
Scores are calculated based on the relative distance of each alternative to the ideal and anti-ideal solutions. Final rankings are determined using the compromise solution, which integrates the positive and negative distances.
}
\value{
A list containing the following components:
  \item{ematrix}{The extended decision matrix, including ideal and anti-ideal solutions.}
  \item{nexmatrix}{The normalized extended decision matrix.}
  \item{wmatrix}{The weighted normalized decision matrix.}
  \item{S}{The aggregated scores for each alternative.}
  \item{k_neg}{The relative closeness to the anti-ideal solution.}
  \item{k_pos}{The relative closeness to the ideal solution.}
  \item{f_k_pos}{The proportion of closeness to the positive ideal solution.}
  \item{f_k_neg}{The proportion of closeness to the negative ideal solution.}
  \item{f_k}{The final compromise score for each alternative.}
  \item{rank}{The ranks of final compromise scores of alternatives.}
}
\references{
Stević, Ž., Pamučar, D., Puška, A., & Chatterjee, P. (2020). Sustainable supplier selection in healthcare industries using a new MCDM method: Measurement of alternatives and ranking according to 
COmpromise solution (MARCOS). \emph{Computers & Industrial Engineering, 140}, 106231. <doş:10.1016/j.cie.2019.106231>
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

# Apply MARCOS with user-defined weights
resmarcos_1 <- marcos(dmatrix=dmat, bcvec=bc, weights=uw)

print(resmarcos_1$f_k)   # Final compromise scores for alternatives
print(resmarcos_1$rank)  # Ranks of the final compromise scores for alternatives
}

