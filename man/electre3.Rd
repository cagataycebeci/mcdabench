\name{electre3}
\alias{electre3}
\encoding{UTF-8}
\title{ELECTRE III (ELimination Et Choix Traduisant la REalité)}
\description{
The ELECTRE III method applies fuzzy outranking by evaluating concordance and discordance between alternatives using user-defined or default pseudo-criteria. It generates a credibility matrix and provides net scores and a partial ranking.
}
\usage{
electre3(dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, 
    thr = 0.5, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric matrix or data frame representing the decision matrix (alternatives as rows, criteria as columns).}
  \item{bcvec}{A numeric vector indicating whether each criterion is a benefit (1) or cost (-1).}
  \item{weights}{A numeric vector of criteria weights, summing to 1 or scaled accordingly.}
  \item{p}{A list of preference thresholds for each criterion. If NULL, default values are estimated.}
  \item{q}{A list of indifference thresholds for each criterion. If NULL, default values are estimated.}
  \item{v}{A list of veto thresholds for each criterion. If NULL, default values are estimated.}
  \item{thr}{The credibility threshold for dominance consideration (default = 0.5).}
  \item{tiesmethod}{A string indicating the method used to resolve ties in ranking. See \code{\link[base]{rank}} for options.}
}
\details{
ELECTRE III uses the pseudo-criteria and fuzzy thresholds.
}
\value{
A list with the following components:
  \item{Credibility}{The n x n credibility matrix among alternatives.}
  \item{Net_Flow}{A numeric vector of average credibility scores per alternative.}
  \item{rank}{A named vector of the ranking of alternatives.}
}
\author{
Cagatay Cebeci
}
\seealso{
\code{\link{electre1}}, \code{\link{electre2}}, \code{\link{electre4}}, \code{\link{calcnormal}}
}
\examples{
dmatrix <- matrix(c(70, 3, 500,
                    90, 4, 400,
                    60, 2, 550), 
                  nrow = 3, byrow = TRUE)
rownames(dmatrix) <- c("A1", "A2", "A3")
bcvec <- c(1, 1, -1)
weights <- c(0.4, 0.3, 0.3)
result <- electre3(dmatrix, bcvec, weights)
result$rank
}
