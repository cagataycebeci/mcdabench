\name{ocra}
\alias{ocra}
\encoding{UTF-8}
\title{OCRA: Operational Competitiveness Rating}
\description{
Computes rankings for decision alternatives using the OCRA method, which evaluates benefit and cost criteria relative to reference values.
}
\usage{
ocra(dmatrix, bcvec, weights, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{A vector indicating the criteria type: 1 for benefit criteria, -1 for cost criteria.}
  \item{weights}{A numeric vector representing the relative importance of each criterion, corresponding to the columns of the decision matrix.}
  \item{tiesmethod}{Method for resolving ties in ranking, defaulting to `"average"`.}
}
\value{
A list with the following components:
  \item{reference_values}{Reference values for each criterion (max for benefit, min for cost).}
  \item{cost_impact_I}{Impact scores for cost criteria (higher values indicate more inefficiency).}
  \item{benefit_impact_O}{Impact scores for benefit criteria (higher values indicate more deviation from ideal benefit levels).}
  \item{total_performance_P}{Final performance scores based on cost and benefit impact values. Lower values indicate better alternatives.}
  \item{rank}{Ranking of alternatives based on performance scores. Lower-ranked alternatives are considered better.}
}
\references{
Madic, M., Petkovic, D., & Radovanovic, M. (2015). Selection of non-conventional machining processes using the OCRA method. \emph{Serbian Journal of Management, 10}(1), 61-73.
}
\examples{
dmatrix <- matrix(c(
  10, 50, 30,  
  20, 60, 20,  
  30, 40, 40,  
  40, 30, 50   
), nrow = 4, byrow = TRUE)

bcvec <- c(1, -1, 1)  # Benefit (1), Cost (-1), Benefit (1)
weights <- c(0.3, 0.4, 0.3)  # Criterion weights

resocra <- ocra(dmatrix, bcvec, weights)
print(resocra)
}
\author{
Cagatay Cebeci
}
