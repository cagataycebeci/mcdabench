\name{aroman}
\alias{aroman}
\encoding{UTF-8}
\title{AROMAN: Alternative Ranking Order Method Accounting for Two Step Normalization}
\description{
  AROMAN applies a two-step normalization process to rank alternatives based on multiple criteria, combining max-min and vector normalization methods.
}
\usage{
aroman(dmatrix, bcvec, weights, beta = 0.5, lambda = 0.5, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numeric matrix (m alternatives x n criteria) representing the initial decision matrix.}
  \item{bcvec}{A numeric vector indicating the type of criteria: \code{1} for benefit criteria (higher is better), \code{-1} for cost criteria (lower is better).}
  \item{weights}{A numeric vector representing the importance of each criterion. If not provided, equal weights are assumed.}
  \item{beta}{A numeric value (default: 0.5) controlling the balance between max-min and vector normalization methods.}
  \item{lambda}{A numeric value (default: 0.5) adjusting the final ranking calculation by weighting cost and benefit criteria.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
 The decision matrix \code{dmatrix} is normalized using both max-min and vector normalization techniques. A weighted normalized matrix is created using the given criteria weights.
 The aggregated average normalization values are calculated.The sum of weighted normalized values is computed separately for benefit and cost criteria.
 A final ranking score (\eqn{R_i}) is calculated using \code{lambda} weighting. Alternatives are ranked based on decreasing \eqn{R_i} values.
}
\value{
  A list containing:
  \item{L}{A numeric vector of weighted sum values for cost criteria.}
  \item{A}{A numeric vector of weighted sum values for benefit criteria.}
  \item{R}{A numeric vector of final ranking values.}
  \item{rank}{An integer vector indicating the ranked order of alternatives.}
}
\references{
Bošković, S., Švadlenka, L., Jovčić, S., Dobrodolac, M., Simić, V., & Bacanin, N. (2023). An alternative ranking order method accounting for two-step normalization (AROMAN) 
- A case study of the electric vehicle selection problem. \emph{IEEE Access, 11}, 39496-39507. <doi:10.1109/ACCESS.2023.3265818>
}
\author{
Cagatay Cebeci
}
\examples{
# Define the decision matrix
dmat <- matrix(c(
  10, 5, 8,
  6, 9, 7,
  8, 7, 9,
  5, 8, 6
), nrow = 4, byrow = TRUE)

# Profit-cost vector
bc <- c(1, -1, 1)

# Weights vector
uw <- c(0.4, 0.3, 0.3)

# Apply AROMAN
resaroman <- aroman(dmatrix = dmat, bcvec = bc, weights = uw, 
    beta = 0.6, lambda = 0.5)
print(resaroman)
}

