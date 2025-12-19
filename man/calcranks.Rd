\name{calcranks}
\alias{calcranks}
\encoding{UTF-8}
\title{Calculate Consensus Ranks using Various Methods}
\usage{
calcranks(rankmatrix, rankmethod = "score")
}
\arguments{
\item{rankmatrix}{A numeric matrix representing the rankings of alternatives by multiple decision-making methods.
Rows correspond to different MCDM methods, and columns correspond to alternatives.
Each cell `rankmatrix[i, j]` should contain the rank of alternative `j` given by method `i` (smaller rank is better).}

\item{rankmethod}{A character string specifying the consensus ranking method to use.
Supported methods are:
\itemize{
  \item \code{"score"}: Calculates the overall score for each alternative based on pairwise comparisons (a higher score indicates a better alternative).
  \item \code{"neustadt"}: Computes Neustadt's rating vector, which is the difference between wins and losses in pairwise comparisons (a higher score indicates a better alternative).
  \item \code{"buchholz"}: Computes Buchholz's rating vector, based on the sum of squared differences of pairwise comparison scores (a lower score indicates a more consistent/less controversial alternative, thus a better rank).
}}
}
\value{
A list containing the results of the ranking calculation:
\itemize{
  \item \strong{\code{rankmethod}}: A character string indicating the ranking method used.
  \item \strong{\code{rating_scores}}: A named numeric vector of the calculated rating scores for each alternative based on the chosen method.
  \item \strong{\code{ordered_alternatives}}: A character vector of alternative names, ordered from the best-performing to the worst-performing according to the calculated \code{rating_scores}.
  \item \strong{\code{score_ranks}}: A numeric vector of the final consensus ranks using the "minimum" ties method.
  \item \strong{\code{average_score_ranks}}: A numeric vector of the final consensus ranks for each alternative the "average" tie method.
}
}
\description{
This function calculates consensus ranks for a set of alternatives, given their rankings from multiple MCDM methods. It supports the Score method, Neustadt's method, and Buchholz's method 
for aggregating ranks.
}
\examples{
# Sample ranking matrix (rows are methods, columns are alternatives R1-R10)
rmat <- matrix(c(
  5, 10, 3, 9, 7, 8, 1, 2, 6, 4,
  10, 9, 6, 3, 7, 1, 5, 2, 4, 8,
  10, 9, 8, 4, 6, 1, 7, 5, 2, 3,
  8, 10, 4, 5, 9, 7, 1, 3, 2, 6,
  9, 7, 1, 8, 3, 4, 2, 6, 10, 5,
  8, 10, 5, 9, 6, 7, 1, 2, 3, 4,
  3, 1, 4, 6, 2, 7, 9, 10, 8, 5,
  2, 1, 5, 3, 4, 6, 8, 10, 9, 7,
  4, 1, 8, 2, 5, 3, 10, 9, 6, 7,
  2, 1, 2, 6, 4, 5, 3, 7, 9, 8,
  2, 1, 7, 4, 3, 5, 10, 9, 8, 6,
  7, 10, 3, 9, 5, 8, 1, 2, 6, 4,
  1, 2, 3, 7, 5, 10, 4, 6, 9, 8,
  4, 1, 6, 3, 2, 5, 10, 9, 8, 7
), nrow = 14, byrow = TRUE)

rownames(rmat) <- c("TOPSIS", "MAUT", "MACBETH", "ORESTE", "VIKOR", "SIR", "EVAMIX", 
   "ARAS", "MOORA", "COPRAS", "WASPAS", "TODIM", "MABAC", "MARE")
colnames(rmat) <- paste0("R", 1:10)

# Calculate ranks using the Score method
results_score <- calcranks(rmat, rankmethod = "score")
print(results_score)

# Calculate ranks using Neustadt's method
results_neustadt <- calcranks(rmat, rankmethod = "neustadt")
print(results_neustadt)

# Calculate ranks using Buchholz's method
results_buchholz <- calcranks(rmat, rankmethod = "buchholz")
print(results_buchholz)
}
\references{
Gogodze, J. (2019). Ranking‐Theory Methods for Solving Multicriteria Decision‐Making Problems. 
 \emph{Advances in Operations Research, 2019}(1), 3217949. <doi:10.1155/2019/3217949>
}
\author{
Cagatay Cebeci
}
\keyword{ranking}
\keyword{mcdm}
\keyword{consensus}
\keyword{decision-making}