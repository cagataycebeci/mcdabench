\name{rankmia}
\alias{rankmia}
\encoding{UTF-8}
\title{Modified Index of Agreement}
\description{
Calculates the Modified Index of Agreement (MIA) between multiple rankings by MCDA methods.
}
\usage{
rankmia(rankmat, j = 1, na.rm = TRUE)
}
\arguments{
\item{rankmat}{A numeric matrix where each row represents a different ranking and each column represents an item being ranked. The values in the matrix are the ranks assigned to each item.}
\item{j}{An integer specifying the exponent for the absolute difference, typically set to 1.}
\item{na.rm}{A logical value indicating whether missing values (\code{NA}) should be removed before calculation. Defaults to \code{TRUE}.}
}
\value{
A data frame where each row and column represents a ranking from the input \code{rankmat}. The values in the data frame are the Modified Index of Agreement scores between the corresponding pairs of rankings, formatted to two decimal places.
}
\details{
The Modified Index of Agreement (MOIA) is an improved version of the Index of Agreement (IA), originally developed by Willmott (1981). It was designed as a standardized measure to evaluate the accuracy of model predictions compared to observed data. 
However, limitations led to modifications that enhance its sensitivity and reliability.

This function extends the application of MOIA to the comparison of rank orderings, useful in Multi-Criteria Decision Analysis (MCDA) problems. The Index of Agreement and Modified Index of Agreement values range between 0 and 1, where 1 indicates 
perfect agreement and 0 indicates complete disagreement.

MIA is calculated using the following formula:

\deqn{MOIA = 1 - \frac{\sum_{i=1}^{N} (|S_i - \bar{O}|^j + |O_i - \bar{O}|^j)}{\sum_{i=1}^{N} |O_i - S_i|^j}}

where:
\describe{
\item{\eqn{MIA}}{The Modified Index of Agreement value (ranging from 0 to 1).}
\item{\eqn{N}}{The total number of items being ranked.}
\item{\eqn{S_i}}{The rank assigned by one MCDA method.}
\item{\eqn{O_i}}{The rank assigned by another MCDA method.}
\item{\eqn{\bar{O}}}{The mean rank value, calculated as \eqn{(N+1)/2}.}
\item{\eqn{j}}{Exponent used for absolute differences (typically \eqn{j=1}).}
\item{\eqn{|\cdot|}}{The absolute value function.}
\item{\eqn{\sum}}{Summation over all ranked items.}
}

Higher values indicate greater agreement between rankings. This formula allows quantification of ranking similarity across MCDA methods.
}
\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},\code{\link{rankentperm}},
  \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}, \code{\link{rankrangesim}}
}

\author{Cagatay Cebeci}
\examples{
# Rank matrix based on MCDA methods
rankmat <- matrix(c(
    3, 12.0, 1, 8, 4.0, 7, 6, 10.0, 2.0, 5, 11.0, 9.0,
    3, 12.0, 2, 8, 4.0, 9, 6, 10.0, 1.0, 5, 11.0, 7.0,
    3, 11.5, 4, 8, 1.5, 7, 5, 9.5, 1.5, 6, 11.5, 9.5,
    4, 3.0, 10, 7, 2.0, 6, 1, 9.0, 5.0, 11, 8.0, 12.0,
    6, 3.0, 10, 5, 2.0, 7, 1, 9.0, 4.0, 11, 8.0, 12.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    10, 1.0, 9, 5, 12.0, 6, 7, 4.0, 11.0, 8, 3.0, 2.0,
    3, 12.0, 1, 8, 4.0, 7, 6, 10.0, 2.0, 5, 11.0, 9.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    5, 3.0, 11, 7, 1.0, 6, 2, 9.0, 4.0, 10, 8.0, 12.0,
    7, 1.0, 11, 9, 2.0, 4, 3, 8.0, 6.0, 10, 5.0, 12.0,
    4, 10.0, 6, 8, 1.0, 5, 3, 9.0, 2.0, 7, 11.0, 12.0,
    3, 9.0, 6, 8, 1.0, 5, 4, 10.0, 2.0, 7, 11.0, 12.0,
    4, 10.0, 6, 8, 1.0, 5, 3, 9.0, 2.0, 7, 11.0, 12.0,
    7, 1.0, 11, 9, 2.0, 4, 3, 8.0, 6.0, 10, 5.0, 12.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    5, 3.0, 11, 7, 2.0, 6, 1, 9.0, 4.0, 10, 8.0, 12.0,
    2, 12.0, 4, 8, 1.0, 7, 6, 9.0, 3.0, 5, 11.0, 10.0,
    1, 12.0, 4, 8, 2.0, 5, 7, 10.0, 3.0, 6, 9.0, 11.0,
    2, 12.0, 1, 8, 4.0, 7, 6, 10.0, 3.0, 5, 11.0, 9.0,
    2, 12.0, 1, 7, 4.0, 8, 6, 10.0, 3.0, 5, 11.0, 9.0,
    4, 7.0, 9, 5, 1.0, 6, 2, 8.0, 3.0, 11, 10.0, 12.0
), nrow=22, byrow=TRUE) # Changed nrow to 22 as there are 22 rows of data

rownames(rankmat) <- c("ARAS", "EDAS", "ELECTRE4", "FUCA", "GRA", "MABAC", "MAIRCA", 
    "MARCOS", "MAUT", "MEGAN", "MOORA", "PROMET1", "PROMET2", "PROMET6", "RAM", "ROV",
    "SMART", "TOPSIS", "VIKOR", "WASPAS", "WPM", "WSM")
colnames(rankmat) <- c("G1", "G2", "G3", "G4", "G5", 
    "G6", "G7", "G8", "G9", "G10", "G11", "G12")

print(rankmat)

# Calculate the Modified Index of Agreement matrix
resmia <- rankmia(rankmat)
print(resmia)
}