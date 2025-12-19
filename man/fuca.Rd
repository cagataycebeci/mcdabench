\name{fuca}
\alias{fuca}
\encoding{UTF-8}
\title{FUCA: Faire Un Choix Adéquat (Make an Adequate Choice)}
\description{
  Implements the Faire Un Choix Adéquat (FUCA) method that ranks alternatives based on the weighted sum of their ranks within each criterion for Multi-Criteria Decision Analysis (MCDA).
}
\usage{
fuca(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows represent alternatives, and columns represent criteria. Row and column names are recommended for better output interpretation.}
  \item{bcvec}{A numeric vector indicating the benefit-cost orientation of each criterion. Use 1 for benefit criteria (higher values are preferred) and -1 for cost criteria (lower values are preferred). 
   The length of this vector must match the number of columns in \code{dmatrix}.}
  \item{weights}{A numeric vector of criteria weights. See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The Faire Un Choix Adéquat (FUCA), which translates from French to "Make an Adequate Choice," is a straightforward and intuitive MCDA method that avoids the need for data normalization (Fernando et al, 2010). 
The algorithm operates by ranking alternatives within each criterion and then aggregating these ranks using criterion weights. 
}
\value{
A list containing the following components:
  \item{Rank_Matrix}{A numeric matrix where each row represents an alternative and each column represents a criterion. The entries are the ranks of the alternatives within each criterion.}
  \item{Weighted_Rank_Scores}{A numeric vector of the overall weighted rank scores for each alternative.}
  \item{Ranking}{A data frame containing the alternatives, their FUCA scores, and their ranks.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower value indicates a better rank).}
}
\references{
Fernando, M. M. L., Escobedo, J. L. P., Azzaro-Pantel, C., Pibouleau, L., Domenech, S., & Aguilar-Lasserre, A. (2011). Selecting the best portfolio alternative from a hybrid multiobjective 
GA-MCDM approach for New Product Development in the pharmaceutical industry. In \emph{2011 IEEE Symposium on Computational Intelligence in Multicriteria Decision-Making (MDCM)} (pp. 159-166). IEEE.
}
\seealso{
  \code{\link{calcweights}}
}
\author{
Cagatay Cebeci
}
\examples{
# Decision matrix (5 alternatives, 4 criteria)
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
userwei <- c(0.1, 0.4, 0.3, 0.2)

# Run FUCA with user-defined weights
resfuca_1 <- fuca(dmat, bcvec=bc, weights=userwei)
print(resfuca_1)

# Run FUCA with the calculated entropy weights
entwei <- calcweights(dmat, bcvec=bc, type="entropy")

resfuca_2 <- fuca(dmat, bcvec=bc, weights=entwei)
print(resfuca_2)
}