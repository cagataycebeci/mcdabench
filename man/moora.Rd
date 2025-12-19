\name{moora}
\alias{moora}
\encoding{UTF-8}
\title{MOORA: Multi-Objective Optimization by Ratio Analysis}
\description{
  The Multi-Objective Optimization by Ratio Analysis (MOORA) method for multi-criteria decision analysis ranks alternatives based on the balance between their performance on benefit and cost criteria.
}
\usage{
moora(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}.}
  \item{weights}{A vector of weights representing the relative importance of each criterion.  See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
  The Multi-Objective Optimization by Ratio Analysis (MOORA) method first normalizes the decision matrix using vector normalization. Then, for each alternative, it calculates an optimization value by summing the weighted normalized values of the benefit criteria and subtracting the weighted normalized values of the cost criteria. Alternatives are ranked based on these optimization values in descending order (higher value implies better performance).
}
\value{
  A list containing the following components:
  \item{Normalized_Matrix}{The normalized decision matrix.}
  \item{Weighted_Normalized_Matrix}{The weighted normalized decision matrix.}
  \item{Optimization_Values}{A numeric vector of the optimization values for each alternative.}
  \item{Ranking}{A data frame with the ranked alternatives, their MOORA scores, and their ranks.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower rank is better).}
}
\references{
Brauers, W. K., & Zavadskas, E. K. (2006). The MOORA method and its application to privatization in a transition economy. \emph{Control and Cybernetics, 35}(2), 445-469.
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

# MOORA with user-defined weights
resmoora_1 <- moora(dmatrix = dmat, bcvec = bc, weights = uw)
print(resmoora_1)

# MOORA with equal weights
equw <- calcweights(dmat, bcvec=bc, type="equal")
resmoora_2 <- moora(dmatrix =  dmat, bcvec = bc, weights = equw)
print(resmoora_2)

# MOORA with entropy weights
entw <- calcweights(dmat, bcvec=bc, type="entropy")
resmoora_3 <- moora(dmatrix = dmat, bcvec = bc, weights = entw)
print(resmoora_3)
}