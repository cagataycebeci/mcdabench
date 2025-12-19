\name{gra}
\alias{gra}
\encoding{UTF-8}
\title{GRA: Grey Relational Analysis}
\description{
  Implements the Grey Relational Analysis (GRA) method for multi-criteria decision analysis by ranking alternatives based on their grey relation degrees with the ideal solution.
}
\usage{
gra(dmatrix, bcvec, weights, idesol = NULL, grdmethod = "sum", 
    rho = 0.5, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{idesol}{An optional numeric vector representing the ideal solution. If \code{NULL} (default), the ideal solution is determined from the normalized matrix (maximum for benefit criteria, minimum for cost criteria). The length must equal the number of columns in \code{dmatrix}.}  
  \item{grdmethod}{The method used to calculate Grey Relational Degrees. Defaults to \code{"sum"}, which applies a weighted sum. If set to \code{"mean"}, a weighted mean is used instead. This choice affects the relative influence of criteria in the final ranking.}
  \item{rho}{A numeric value between 0 and 1 (default is 0.5) representing the distinguishing coefficient. It influences the grey relational coefficients.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
  The Grey Relational Analysis (GRA) method normalizes the decision matrix and then calculates the grey relational coefficients between each alternative and the ideal solution. The grey relation degree for each alternative is the average of its grey relational coefficients across all criteria. Alternatives are ranked based on their grey relation degrees in descending order (higher degree implies closer to the ideal solution). The grey relational coefficients for each criterion are multiplied by their respective weights before calculating the grey relation degree. To calculate Grey Relational Degrees, \code{"sum"} applies a weighted sum. If \code{grdmethod} is set to \code{"mean"}, a weighted mean is used instead. This choice affects the relative influence of criteria in the final ranking.
}
\value{
  A list containing the following components:
  \item{Ideal_Solution}{The ideal solution vector used for comparison.}
  \item{Grey_Relation_Coef}{A numeric matrix of grey relational coefficients for each alternative and criterion.}
  \item{Grey_Relation_Deg}{A numeric vector of grey relation degrees for each alternative.}
  \item{Ranking}{A data frame with the ranked alternatives, their grey relation degrees, and their ranks.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower rank is better).}
}
\seealso{
    \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\references{
Julong, D. (1989). Introduction to grey system theory. \emph{The Journal of Grey System, 1}(1), 1-24.
Deng, J.L. (1982), Control problems of grey systems, \emph{Systems and Control Letters, 5}, 288‐294.
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

# Apply GRA with default settings
resgra_1 <- gra(dmatrix=dmat, bcvec=bc, weights=uw)
print(resgra_1)

# Apply GRA with a low rho value
resgra_2 <- gra(dmat, bc, weights = uw, rho = 0.3)
print(resgra_2)

# Apply GRA with a higher rho value
resgra_3 <- gra(dmatrix=dmat, bcvec=bc, weights=uw, grdmethod="sum", rho = 0.7)
print(resgra_3)

# Apply GRA with ideal solutions
ideal_solutions <- c(80, 15, 6, 16)
resgra_4 <- gra(dmatrix=dmat, bcvec=bc, weights=uw, idesol=ideal_solutions, rho=0.5)
print(resgra_4)
}