\name{cocoso}
\alias{cocoso}
\encoding{UTF-8}
\title{COCOSO: Combined Compromise Solution Method}
\description{
  COCOSO proposed by Yazdani et al (2019) is an MCDA method that normalizes the decision matrix and ranks alternatives based on specified criterion weights.
}
\usage{
  cocoso(dmatrix, bcvec, weights, lambda = 0.5, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A \code{matrix} containing the decision matrix. Alternatives are in rows, criteria are in columns.}
  \item{bcvec}{A vector indicating the preferred direction for each criterion (\code{c(1, -1, ...)}). Positive values (\code{1}) mean higher values are preferred, while negative values (\code{-1}) mean lower values are preferred.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{lambda}{A weighting parameter (\code{0.5} by default).}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\value{
  Returns a list containing:
  \item{S_i}{Total weighted normalized values of alternatives.}
  \item{P_i}{Multiplicative values of alternatives.}
  \item{k_ia}{First combined calculated value.}
  \item{k_ib}{Second combined calculated value.}
  \item{k_ic}{Third combined calculated value.}
  \item{k_i}{Final composite value.}
  \item{rank}{Ranking of the alternatives.}
}
\author{
  Cagatay Cebeci
}
\references{
Yazdani, M., Zarate, P., Kazimieras Zavadskas, E., & Turskis, Z. (2019). A combined compromise solution (CoCoSo) method for multi-criteria decision-making problems.
\emph{Management Decision, 57}(9), 2501-2519. <doi:10.1108/MD-05-2017-0458>
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
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

# Apply COCOSO with user-defined weights
rescocoso_1 <- cocoso(dmatrix=dmat, bcvec = bc, weights = uw, lambda=0.1)
print(rescocoso_1)

# Apply COCOSO with CRITIC weights
cw <- calcweights(dmatrix=dmat, bcvec=bc, type="critic", normethod="linear")
rescocoso_2 <- cocoso(dmatrix=dmat, bcvec = bc, weights = cw, lambda=0.8, tiesmethod="average")
print(rescocoso_2)
}