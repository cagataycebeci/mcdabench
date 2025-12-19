\name{macont6}
\alias{macont6}
\encoding{UTF-8}
\title{MACONT T6: Multi-Criteria Decision-Making with Comprehensive Normalization (6 steps)}
\description{
MACONT 6 applies a multi-step normalization process to rank alternatives based on multiple criteria using a combination of arithmetic and geometric weighted averaging.
}
\usage{
macont6(dmatrix, bcvec, weights, p=0.5, q=0.5, delta=0.5, theta=0.5, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix (m alternatives x n criteria) representing the initial decision matrix.}
  \item{bcvec}{A numeric vector indicating the type of criteria: \code{1} for benefit criteria (higher is better), \code{-1} for cost criteria (lower is better).}
  \item{weights}{A numeric vector representing the importance of each criterion. If not provided, equal weights are assumed.}
  \item{p}{A numeric value (between 0 and 1) representing the weight for the arithmetic averaging operator.}
  \item{q}{A numeric value (between 0 and 1) representing the weight for the geometric averaging operator.}
  \item{delta}{A numeric value (between 0 and 1) controlling the influence of rho and zeta in the ranking process.}
  \item{theta}{A numeric value (between 0 and 1) adjusting preference functions for ranking.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
  \enumerate{
    \item The decision matrix \code{dmatrix} is normalized using three different techniques: ratio-based, max-min, and vector normalization.
    \item A balanced normalized matrix is computed using weighted averages of the three normalization methods.
    \item The S1 and S2 ranking values are calculated using rho and zeta components.
    \item A final ranking score (\eqn{S}) is derived, adjusting for preference weighting.
    \item Alternatives are ranked in descending order based on comprehensive scores.
  }
}
\value{
  A list containing:
  \item{balance_matrix}{A numeric matrix of balanced normalization values.}
  \item{S1_values}{A numeric vector of S1 ranking values.}
  \item{S2_values}{A numeric vector of S2 ranking values.}
  \item{comprehensive_scores}{A numeric vector of final ranking values (\eqn{S}).}
  \item{rank}{An integer vector indicating the ranked order of alternatives.}
}
\references{
Wen, Z., Liao, H., & Zavadskas, E. K. (2020). MACONT: Mixed aggregation by comprehensive normalization technique for multi-criteria analysis. \emph{Informatica, 31}(4), 857-880. <doi:10.15388/20-INFOR417>

Nguyen, A. T. (2023). Expanding the data normalization strategy to the MACONT method for multi-criteria decision making. \emph{Engineering, Technology & Applied Science Research, 13}(2), 10489-10495. <doi:10.48084/etasr.5672>
}
\seealso{
 \code{\link{topsis}}, \code{\link{vikor}}, \code{\link{electre1}}
}
\author{
Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(
  5, 3, 4,
  4, 5, 2,
  3, 4, 5,
  5, 2, 3
), nrow = 4, byrow = TRUE)

bcvec <- c(1, -1, 1)
weights <- c(0.3, 0.4, 0.3)

# Run MACONT
res_macont <- macont6(dmatrix = dmat, bcvec = bcvec, weights = weights,
     p = 0.4, q = 0.3, delta = 0.6, theta = 0.7)

print(res_macont$rank)
}

