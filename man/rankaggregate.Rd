\name{rankaggregate}
\alias{rankaggregate}
\encoding{UTF-8}
\title{Aggregat Ranks of Alternatives for Decision Making}
\description{
  Aggregates and summarizes the ranking information of alternatives obtained from 
  multiple methods. Given a rank matrix where rows represent evaluation methods 
  (e.g., MCDA methods) and columns represent alternatives, the function computes 
  aggregated overall ranks, produces hierarchical outranking flow, 
  and counts the frequency of each alternative appearing in the top \code{tk} ranks.
}
\usage{
rankaggregate(rankmat, topk=3, damp=0.5, niters=100, tol=1e-4, tiesmethod="average")
}
\arguments{
  \item{rankmat}{A numeric matrix in which rows correspond to different ranking 
    methods and columns correspond to alternatives. Each element \code{rankmat[i, j]} 
    indicates the rank of alternative \code{j} based on method \code{i}, with lower 
    values indicating better performance.}
  \item{topk}{A positive integer indicating the "top \code{tk}" positions. For each method, 
    the function counts how many times each alternative appears within the top \code{tk} ranks. 
    The default is 3.}
  \item{damp}{A numeric value in \eqn{(0, 1)} used as the damping factor in the Markov 
    chain aggregation method. The default is 0.5.}
  \item{niters}{A positive integer specifying the maximum number of iterations allowed in 
    the Markov chain aggregation method. The default is 200.}
  \item{tol}{A numeric value specifying the convergence tolerance for the iterative 
    procedure used in the Markov chain aggregation method. The default is 1e-4.}
  \item{tiesmethod}{A character string specifying the method used to handle ties when 
    computing overall ranks. See \code{\link[base]{rank}} for available options 
    (e.g., \code{"average"}, \code{"first"}, \code{"random"}, \code{"max"}, \code{"min"}). 
    The default is \code{"average"}.}
}
\details{
  The \code{prefsummary} function aggregates ranking information for a set of alternatives 
  evaluated by multiple methods. In particular, it performs the following tasks:
  \enumerate{
    \item Computes alternative rankings based on the Top-\code{tk} counts method. 
    \item Computes the overall ranks based on the sum of the ranks across methods using the Rank Sum approach.
    \item Produces aggregated ranks using the Median of ranks method.
    \item Computes pairwise comparisons between alternatives based on the Copeland method 
      (assigning 1 point for a win and 0.5 points for a tie in pairwise comparisons).
    \item Calculates overall scores using the Kemeny & Young method based on Kendall Tau distances.
    \item Computes a stationary score and ranking using a Markov chain (PageRank-like) 
      approach.
    \item Combines the rankings from the above methods into a matrix and constructs a 
      hierarchical (outranking) string (preference hierarchy) for each method.
  }
}
\value{
  A list with the following components:
  \item{restopk}{A list from the Top-\code{tk} counts method containing the top-\code{tk} 
    threshold, counts for each alternative, and the corresponding ranking.}
  \item{resranksum}{A list containing the rank sums for each alternative and the ranking based on these sums.}
  \item{resmedian}{A list containing the median rank for each alternative and its ranking.}
  \item{rescopeland}{A list containing the Copeland scores for each alternative and the corresponding ranking.}
  \item{reskemeny}{A list containing the Kemeny & Young scores for each alternative and the corresponding ranking.}
  \item{resmarkov}{A list containing the Markov chain results: number of iterations, the wins matrix, the transition matrix, the Markov scores, and the resulting ranking.}
  \item{preference_ranking}{A matrix that aggregates the rankings from all methods. Row names correspond to methods while column names correspond to alternatives.}
  \item{preference_table}{A data frame presenting, for each method, a hierarchical outranking string (i.e., a preference hierarchy among alternatives).}
}
\references{
Saari, D. G., & Merlin, V. R. (1996). The Copeland method: I.: Relationships and the dictionary. \emph{Economic Theory, 8}, 51-76. <doi:10.1007/BF01212012>
}
\author{
Cagatay Cebeci
}
\examples{
# Example rank matrix (rows: methods, columns: alternatives)
ranks <- c(
  5.0, 7.0, 6.0, 10.0, 8.0, 9.0, 3.0, 1.0, 4.0, 2.0,
  7.0, 8.0, 3.5, 3.5, 9.0, 3.5, 3.5, 3.5, 3.5, 10.0,
  5.0, 7.0, 6.0, 10.0, 9.0, 8.0, 3.0, 2.0, 4.0, 1.0,
  6.0, 5.0, 3.0, 9.0, 8.0, 10.0, 1.0, 2.0, 7.0, 4.0,
  5.5, 7.0, 4.0, 10.0, 9.0, 8.0, 3.0, 1.0, 5.5, 2.0,
  9.0, 7.0, 6.0, 8.0, 2.0, 10.0, 4.0, 5.0, 1.0, 3.0,
  10.0, 9.0, 6.0, 7.0, 4.0, 8.0, 5.0, 3.0, 2.0, 1.0,
  5.0, 7.0, 6.0, 10.0, 8.0, 9.0, 3.0, 1.0, 4.0, 2.0,
  5.0, 7.0, 6.0, 10.0, 8.0, 9.0, 3.0, 1.0, 4.0, 2.0,
  5.0, 7.0, 6.0, 10.0, 8.0, 9.0, 3.0, 1.0, 4.0, 2.0,
  2.0, 8.0, 3.0, 9.0, 5.0, 10.0, 1.0, 6.0, 4.0, 7.0,
  2.0, 8.0, 3.0, 9.0, 5.0, 10.0, 1.0, 6.0, 4.0, 7.0,
  9.0, 6.0, 7.0, 5.0, 3.0, 8.0, 4.0, 2.0, 2.0, 1.0,
  10.0, 6.0, 8.0, 7.0, 5.0, 9.0, 4.0, 3.0, 2.0, 1.0,
  10.0, 6.0, 7.0, 9.0, 5.0, 8.0, 4.0, 3.0, 2.0, 1.0,
  7.0, 10.0, 5.0, 9.0, 6.0, 8.0, 3.0, 1.0, 4.0, 2.0,
  6.0, 9.0, 5.0, 10.0, 7.0, 8.0, 3.0, 1.0, 4.0, 2.0,
  7.0, 10.0, 5.0, 9.0, 6.0, 8.0, 3.0, 1.0, 4.0, 2.0,
  3.0, 3.0, 10.0, 3.0, 3.0, 3.0, 8.0, 6.0, 9.0, 7.0,
  3.0, 9.0, 5.0, 10.0, 8.0, 7.0, 6.0, 2.0, 4.0, 1.0,
  10.0, 6.0, 8.0, 7.0, 5.0, 9.0, 4.0, 3.0, 2.0, 1.0,
  2.0, 8.0, 3.0, 9.0, 5.0, 10.0, 1.0, 6.0, 4.0, 7.0,
  2.0, 8.0, 3.0, 9.0, 5.0, 10.0, 1.0, 6.0, 4.0, 7.0,
  8.0, 6.0, 5.0, 10.0, 7.0, 9.0, 3.0, 1.0, 4.0, 2.0,
  9.0, 3.0, 7.0, 5.0, 10.0, 1.0, 8.0, 2.0, 6.0, 4.0,
  8.0, 9.0, 3.0, 4.0, 10.0, 6.0, 2.0, 1.0, 5.0, 7.0,
  8.5, 8.5, 3.0, 5.0, 8.5, 6.0, 2.0, 1.0, 4.0, 8.5,
  5.0, 7.0, 6.0, 10.0, 8.0, 9.0, 3.0, 1.0, 4.0, 2.0
)

rankmat <- matrix(ranks, nrow = 28, byrow = TRUE)
rownames(rankmat) <- c(
  "ARAS", "COCOSO", "CODAS", "EDAS", "ELECTRE4", "FUCA", "GRA", "MABAC",
  "MAIRCA", "MARCOS", "MAUT", "MAVT", "MEGAN", "MOORA", "PROMET1",
  "PROMET2", "PROMET3", "PROMET4", "PROMET5", "PROMET6", "RAM", "ROV",
  "SMART", "TOPSIS", "VIKOR", "WASPAS", "WPM", "WSM"
)
colnames(rankmat) <- paste0("ALT_", 1:ncol(rankmat))

# Aggregate the ranks
res_summary <- rankaggregate(rankmat, tiesmethod = "average", 
   topk = 1, tol = 1e-4, niters = 10, damp = 0.4)
print(res_summary)
}
