\name{electre2}
\alias{electre2}
\encoding{UTF-8}
\title{ELECTRE II (ELimination Et Choix Traduisant la REalité)}
\description{
Implements the ELECTRE II outranking method to rank a set of alternatives based on concordance, discordance, and credibility analysis.
}
\usage{
electre2(dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, 
         st = 0.65, wt = 0.5, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric matrix or data frame representing the decision matrix, where rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the nature of each criterion: \code{1} for benefit and \code{-1} for cost.}
  \item{weights}{A numeric vector containing weights for each criterion. Must match the number of columns in \code{dmatrix}.}
  \item{p}{A list of preference thresholds for each criterion. If \code{NULL}, default values are estimated from the data.}
  \item{q}{A list of indifference thresholds for each criterion. If \code{NULL}, defaults are estimated.}
  \item{v}{A list of veto thresholds for each criterion. If \code{NULL}, defaults are estimated.}
  \item{st}{Strong credibility threshold (default = 0.65). Values above this define a strong outranking relation.}
  \item{wt}{Weak credibility threshold (default = 0.5). Values between \code{wt} and \code{st} define a weak outranking relation.}
  \item{tiesmethod}{String specifying the tie-breaking method to be used in \code{\link[base]{rank}}.}
}
\details{
ELECTRE II uses pairwise credibility analysis with two credibility thresholds to distinguish strong and weak outranking relations among alternatives.

The \code{st} and \code{wt} arguments define the strong and weak credibility cutoffs. While there is no fixed default specified in the original ELECTRE II methodology, values such as \code{st = 0.65} and \code{wt = 0.5} are commonly used in the literature due to their balance between discrimination and robustness.

However, users are encouraged to tune these thresholds based on the decision context, the number of alternatives, and how sensitive the comparison process should be. In highly uncertain decision environments, lower thresholds may be appropriate.

The final ranking is computed as the average rank between upward and downward outranking scores, providing a balanced assessment of dominance and vulnerability.
}
\value{
A list with the following elements:
  \item{Concordance}{Matrix of concordance indices between alternatives.}
  \item{Discordance}{Matrix of discordance indices between alternatives.}
  \item{Credibility}{Overall credibility matrix based on concordance and discordance.}
  \item{Strong}{Logical matrix indicating strong outranking relationships.}
  \item{Weak}{Logical matrix indicating weak outranking relationships.}
  \item{Rank_Upward}{Ranking of alternatives from strongest outranking profile.}
  \item{Rank_Downward}{Ranking of alternatives from most dominated profile.}
  \item{rank}{Final overall ranking based on average of upward and downward ranks.}
}
\examples{
dmatrix <- matrix(c(70, 3, 500,
                    90, 4, 400,
                    60, 2, 550),
                  nrow = 3, byrow = TRUE)
rownames(dmatrix) <- c("A1", "A2", "A3")
bcvec <- c(1, 1, -1)
weights <- c(0.4, 0.3, 0.3)
result <- electre2(dmatrix, bcvec, weights)
result$rank
}
\seealso{
\code{\link{electre3}}, \code{\link{electre4}}, \code{\link{calcnormal}}, \code{\link[base]{rank}}
}
\author{
Cagatay Cebeci
}
