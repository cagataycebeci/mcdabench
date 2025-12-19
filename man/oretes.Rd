\name{oretes}
\alias{oretes}
\encoding{UTF-8}
\title{ORETES: Organísation, Rangemént et Estímatíon En Synthèse}
\description{
  The "Organísation, Rangemént et Estímatíon En Synthèse (ORETES)" is an MCDM technique based on outranking relations and graph theory.The function calculates lambda values, psi ranks, a normalized dominance matrix, outranking flows, and provides a final ranking of alternatives. It also offers an optional visualization of the dominance graph.
}
\usage{
oretes(dmatrix, bcvec, weights, tiesmethod = "average", domplot = FALSE)
}
\arguments{
  \item{dmatrix}{A numeric matrix or data frame representing the decision matrix.
    Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector of the same length as the number of criteria,
    indicating the type of each criterion. Use \code{1} for benefit criteria (higher is better)
    and \code{-1} for cost criteria (lower is better).}
  \item{weights}{A numeric vector of the same length as the number of criteria,
    representing the weights (importance) of each criterion. These do not necessarily need to sum to 1.}
  \item{tiesmethod}{A character string specifying how ties are handled when ranking.
    Passed to \code{\link[base]{rank}}. Common options include \code{"average"} (default),
    \code{"first"}, \code{"last"}, \code{"random"}, \code{"max"}, \code{"min"}.}
  \item{domplot}{A logical value. If \code{TRUE}, the dominance graph illustrating
    the outranking relationships between alternatives will be plotted. Requires the \code{igraph} package.}
}
\details{
  The ORETES method proceeds through several steps:
  \enumerate{
    \item \bold{Lambda Values Calculation}: Normalizes criterion values for each alternative
          to a scale of 0 to 1, where 1 represents the ideal performance and 0 the worst.
    \item \bold{Psi Matrix Calculation}: Ranks the lambda values for each alternative across
          its criteria. This provides an internal ranking of criteria importance for each alternative.
    \item \bold{Dominance Matrix Calculation}: Computes the strength of preference (dominance)
          of one alternative over another. For each pair of alternatives (A, B), it sums the weights
          of the criteria where A performs better than B. The matrix is then normalized.
    \item \bold{Outranking Flow Calculation}: Determines the net outranking flow for each alternative
          by subtracting its incoming dominance strength from its outgoing dominance strength.
          A higher net flow indicates a stronger overall preference.
    \item \bold{Final Ranking}: Alternatives are ranked based on their net outranking flow,
          from highest to lowest.
  }
  The dominance graph, if plotted, visually represents the normalized dominance matrix,
  with arrows indicating the direction and strength of preference between alternatives.
  Thicker arrows indicate stronger dominance.
}
\value{
  A list containing the following components:
  \item{lambda_matrix}{A numeric matrix of normalized lambda values (0-1) for each alternative across each criterion.}
  \item{psi_matrix}{A numeric matrix where each row represents an alternative and contains the ranks of its lambda values across criteria.}
  \item{dominance_matrix}{A numeric matrix showing the normalized strength of dominance of row alternatives over column alternatives.}
  \item{outranking_flow}{A numeric vector containing the net outranking flow for each alternative.}
  \item{rank}{A numeric vector indicating the final rank of each alternative, where 1 is the best rank.}
}
\references{
Roubens, M. (1982). Preference relations on actions and criteria in multicriteria decision making. European Journal of Operational Research, 10(1), 51-55. <doi:10.1016/0377-2217(82)90131-X>
}
\author{
  Cagatay Cebeci
}
\note{
  This implementation of ORETES is an adaptation focused on the core principles of
  lambda values, dominance relations, and outranking flow. It simplifies some of
  the more intricate aspects of the original method, such as specific indifference/preference
  thresholds and detailed preference functions, for broader applicability in MCDM.
}
\seealso{
  \code{\link[igraph]{graph_from_data_frame}}, \code{\link[igraph]{plot.igraph}}
}
\examples{
# Sample decision matrix
dmatrix <- matrix(c(
  10, 50, 30,  
  20, 60, 20, 
  30, 40, 40, 
  40, 30, 50), 
  nrow = 4, byrow = TRUE, 
  dimnames = list(paste0("A", 1:4), paste0("C", 1:3)
))

# Define benefit (1) or cost (-1) for each criterion
bcvec <- c(1, -1, 1) 

# Define weights for each criterion
weights <- c(0.3, 0.4, 0.3) # Sum does not necessarily need to be 1

# ORETES without plotting the dominance graph
resoreste_1 <- oretes(dmatrix, bcvec, weights, domplot = FALSE)
print(resoreste_1)

# ORETES with the dominance graph
resoreste_2 <- oretes(dmatrix, bcvec, weights, domplot = TRUE)
print(resoreste_2)

# Zero weights (no dominance relationships)
  zerowei <- c(0, 0, 0)
  resoreste_3 <- oretes(dmatrix, bcvec, zerowei, domplot = TRUE)
  print(resoreste_3)
}
