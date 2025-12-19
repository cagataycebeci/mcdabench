\name{edas}
\alias{edas}
\title{EDAS: Evaluation based on Distance from Average Solution}
\description{
This function implements the method Evaluation based on Distance from Average Solution (EDAS) (Keshavarz et al, 2015) for multi-criteria decision-making by calculating distances from the average solution for each criterion, and evaluating alternatives based on positive and negative distances.
}
\usage{
edas(dmatrix, bcvec, weights, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{Decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).}
  \item{weights}{A numeric vector of criteria weights. See \code{\link{calcweights}} for details.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The EDAS method evaluates alternatives based on their positive and negative distances from the average solution for each criterion. Positive distances (\code{pda}) measure how much better an alternative is compared to the average, while negative distances (\code{nda}) measure how much worse it is. These values are combined into scores for ranking alternatives.
}
\value{
A list containing the following components:
  \item{amatrix}{The average solution for each criterion.}
  \item{pda}{The positive distance from the average for each alternative and criterion.}
  \item{nda}{The negative distance from the average for each alternative and criterion.}
  \item{sp}{The aggregated positive distances for each alternative.}
  \item{sn}{The aggregated negative distances for each alternative.}
  \item{nsp}{The normalized positive distances for each alternative.}
  \item{nsn}{The normalized negative distances for each alternative.}
  \item{score}{The final scores of alternatives, based on positive and negative distances.}
  \item{rank}{The ranks of final score of alternatives.}
}
\references{
Keshavarz Ghorabaee, M., Zavadskas, E. K., Olfat, L., & Turskis, Z. (2015). Multi-criteria inventory classification using a new method of evaluation based on distance from average solution (EDAS). \emph{Informatica, 26(3)}, 435-451.<doi: 10.3233/INF-2015-1070>.
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

# Apply EDAS with user-defined weights
resedas_1 <- edas(dmatrix=dmat, bcvec = bc, weights = userwei)
print(resedas_1)

# Apply EDAS with GINI weights
giniwei <- calcweights(dmatrix=dmat, bcvec = bc, type="gini")
resedas_2 <- edas(dmatrix=dmat, bcvec = bc, weights = giniwei)
print(resedas_2) 
}
