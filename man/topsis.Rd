\name{topsis}
\alias{topsis}
\encoding{UTF-8}
\title{TOPSIS: Technique for Order of Preference by Similarity to Ideal Solution}
\description{
This function implements the Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) method for multi-criteria decision-making by evaluating alternatives based on their distances from the ideal and anti-ideal solutions.
}
\usage{
topsis(dmatrix, bcvec, weights, normethod = "maxmin", tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{
    Decision matrix where rows represent alternatives and columns represent criteria.
  }
  \item{bcvec}{
    Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).
  }
  \item{weights}{A numeric vector of user-defined weights representing the relative importance of each criterion, or a character string to calculate it internally. Defaults to \code{"equal"} to apply equal weighting for criteria. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"maxmin"}. See \code{\link{calcnormal}} for details. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) method for Multi-Criteria Decision Making (MCDM) was originally developed by Ching-Lai Hwang and Kwangsun Yoon in 1981. Their seminal work, "Multiple Attribute Decision Making: Methods and Applications", published by Springer-Verlag in 1981, introduced the fundamental concepts and methodology of TOPSIS. Later, there were further developments and refinements of the method by Yoon in 1987, and by Hwang, Lai, and Liu in 1993.
The TOPSIS method evaluates alternatives by normalizing the decision matrix, applying weights, and calculating distances to the positive ideal solution (\code{pis}) and the negative ideal solution (\code{nis}). The closeness coefficient (\code{p}) is computed as the ratio of the distance to the negative ideal solution (\code{Dm}) to the sum of the distances to both ideal solutions (\code{Dm + Dp}).
}
\value{
A list containing the following components:
  \item{normalized_matrix}{The normalized decision matrix.}
  \item{weighted_matrix}{The weighted decision matrix.}
  \item{pis}{The positive ideal solution for each criterion.}
  \item{nis}{The negative ideal solution for each criterion.}
  \item{Dm}{The distance from each alternative to the negative ideal solution.}
  \item{Dp}{The distance from each alternative to the positive ideal solution.}
  \item{p}{The closeness coefficient for each alternative.}
  \item{rank}{The rank of closeness coefficients.}
}
\references{
Hwang, C.L. & Yoon, K. (1981). Multiple Attribute Decision Making: Methods and Applications. New York: Springer-Verlag.

Yoon, K. (1987). A reconciliation among discrete compromise solutions. \emph{Journal of the Operational Research Society, 38}(3), 277-286. <doi:10.1057/jors.1987.44>

Hwang, C. L., Lai, Y. J., & Liu, T. Y. (1993). A new approach for multiple objective decision making. \emph{Computers & Operations Research, 20}(8), 889-899.

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

# Apply TOPSIS with user-defined weights
restopsis_1 <- topsis(dmatrix = dmat, bcvec = bc, weights = uw)
print(restopsis_1) 

# Apply TOPSIS with entropy weighting method
entw <- calcweights(dmat, bcvec=bc, type="entropy")
restopsis_2 <- topsis(dmatrix = dmat, bcvec = bc, weights = entw)
print(restopsis_2$p) # Closeness coefficients for ranking alternatives
print(restopsis_2$rank) # Rank of alternatives
}
