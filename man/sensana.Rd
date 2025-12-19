\name{sensana}
\alias{sensana}
\encoding{UTF-8}
\title{Sensitivity and Stability Analyses for MCDA}
\description{
Performs various sensitivity and stability analyses on a matrix of rank results from Multi-Criteria Decision Analysis (MCDA).
}
\usage{
sensana(rankmat)
}
\arguments{
  \item{rankmat}{A numeric matrix where rows represent alternatives and columns represent different scenarios (e.g., different weighting vectors or MCDA methods). Each element in the matrix is the rank of the alternative in the corresponding scenario.}
}
\details{
This function calculates several measures to assess the sensitivity and stability of alternative rankings across different scenarios. The computed measures include:

\itemize{
  \item{\code{SD}:} The standard deviation of the ranks for each alternative across all scenarios, indicating the dispersion of its ranks.
  \item{\code{RSI}:} Rank Stability Index (Stability scores) for each alternative, calculated as 1 minus the proportion of rank changes across scenarios. A higher score indicates greater rank stability.
  \item{\code{RVOL}:} Rank Volatility scores for each alternative, representing the proportion of rank changes between consecutive scenarios. A lower score indicates less rank volatility.
  \item{\code{sensscores}:} Weight Sensitivity scores for each scenario, representing the average absolute rank change of all alternatives when compared to the first scenario. Lower scores suggest higher stability of the overall ranking relative to the first scenario.
}
}
\value{
A list containing the following components:
  \item{stability_table}{A data frame of stability measures for each alternative.}
  \item{spearmancor}{A numeric matrix of Spearman's rank correlation coefficients between scenarios (between -1 and 1).}
  \item{sensscores}{A numeric vector of weight sensitivity scores for each scenario.}
}

\author{
Cagatay Cebeci
}
\examples{
# Example rank matrix
rankmat <- matrix(
  c(1, 2, 1, 3,
    2, 1, 3, 1,
    3, 3, 2, 2),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("ALT1", "ALT2", "ALT3"), paste0("S_", 1:4))
)

# Perform sensitivity and stability analysis
results <- sensana(rankmat)

print(results)
}