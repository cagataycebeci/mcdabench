\name{codas}
\alias{codas}
\encoding{UTF-8}
\title{CODAS: COmbinative Distance-based ASsessment}
\description{
  Implements the COmbinative Distance-based ASsessment (CODAS) method for multi-criteria decision analysis (Keshavarz et al, 2016).
  This method ranks alternatives based on a combination of Euclidean and Taxicab distances from the negative ideal solution.
}
\usage{
codas(dmatrix, bcvec, weights, thr = 0.1, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{thr}{A numeric value between 0 and 1 representing the threshold parameter. It influences the balance between Euclidean and Taxicab distances in the appraisal score calculation. Default is 0.1.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
  The CODAS method evaluates alternatives by calculating their distances from the negative ideal solution using both Euclidean and Taxicab metrics. A combinative appraisal score is then computed, 
  giving more preference to the Euclidean distance. The alternatives are ranked based on these appraisal scores in descending order.
}
\value{
  A list containing the following components:
  \item{niv}{The negative ideal solution vector.}
  \item{euc_dist}{A numeric vector of Euclidean distances of each alternative from the negative ideal solution.}
  \item{tax_dist}{A numeric vector of Taxicab distances of each alternative from the negative ideal solution.}
  \item{pi_i}{A numeric vector of appraisal scores for each alternative.}
  \item{rank}{A numeric vector representing the rank of each alternative. Lower rank is better.}
}

\references{
  Keshavarz Ghorabaee, M., Zavadskas, E. K., Turskis, Z., & Antuchevičienė, J. (2016). A new combinative distance-based assessment (CODAS) 
  method for multi-criteria decision-making. <https://etalpykla.vilniustech.lt/handle/123456789/116529>.
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
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
uw <- c(0.1, 0.4, 0.3, 0.2)

# Apply CODAS with user-defined weights, maxmin normalization, and thr = 0.1
rescodas_1 <- codas(dmatrix=dmat, bcvec = bc, weights = uw)
print(rescodas_1)

# Apply CODAS with entropy weights
entwei <- calcweights(dmatrix=dmat, bcvec = bc, type="entropy", normethod="sum")
rescodas_2 <- codas(dmatrix=dmat, bcvec = bc, weights = entwei, thr = 0.5)
print(rescodas_2)
}