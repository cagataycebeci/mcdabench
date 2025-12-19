\name{electre1}
\alias{electre1}
\encoding{UTF-8}
\title{ELECTRE I (ELimination Et Choix Traduisant la REalité)}
\description{
Performs outranking analysis using the ELECTRE I method, identifying dominant alternatives based on concordance and discordance indices.
}
\usage{
electre1(dmatrix, bcvec, weights, ct = 0.65, dt = 0.35, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{A numeric vector indicating for each criterion whether higher values are preferable (\code{1}) or not (\code{-1}).}
  \item{weights}{A numeric vector of weights assigned to each criterion. Length must match the number of columns in \code{dmatrix}.}
  \item{ct}{Concordance threshold. Must be between 0 and 1. Defaults to \code{0.65}.}
  \item{dt}{Discordance threshold. Must be between 0 and 1. Defaults to \code{0.35}.}
  \item{tiesmethod}{Tie-breaking method for final ranking. Passed to \code{\link[base]{rank}}. Defaults to \code{"average"}.}
}
\details{
ELECTRE I is the original outranking-based multi-criteria decision method. It constructs binary pairwise relations based on:

\itemize{
  \item Concordance index: weighted agreement over criteria supporting the preference of one alternative over another.
  \item Discordance index: the strongest opposition to that preference.
}
An alternative \(a_i\) outranks \(a_j\) only if:

\itemize{
  \item The concordance index is at least \code{ct}, and
  \item The discordance index is at most \code{dt}.
}
The final result includes the outranking relations and the kernel (alternatives that are not strongly dominated by others). 
Additionally, a basic ranking is computed based on how many alternatives each option outranks, using \code{tiesmethod}.
}

\value{
A list with the following elements:
  \item{Concordance}{The n x n matrix of concordance indices.}
  \item{Discordance}{The n x n matrix of discordance indices.}
  \item{Outrank}{Logical matrix indicating outranking relations between alternatives.}
  \item{Kernel}{Indices of non-dominated alternatives (the kernel).}
  \item{Outrank_Score}{Number of alternatives outranked by each alternative.}
  \item{rank}{Ranking based on outrank strength.}
}
\examples{
dmatrix <- matrix(c(70, 3, 500,
                    90, 4, 400,
                    60, 2, 550), 
                  nrow = 3, byrow = TRUE)
rownames(dmatrix) <- c("A1", "A2", "A3")
bcvec <- c(1, 1, -1)
weights <- c(0.4, 0.3, 0.3)
result <- electre1(dmatrix, bcvec, weights)
result$Kernel
result$rank
}
\seealso{
\code{\link{electre2}}, \code{\link{electre3}}, \code{\link{electre4}}, \code{\link{calcnormal}}
}
\author{
Cagatay Cebeci
}
