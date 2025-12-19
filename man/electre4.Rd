\name{electre4}
\alias{electre4}
\encoding{UTF-8}
\title{ELECTRE IV (ELimination Et Choix Traduisant la REalité)}
\description{
The ELECTRE IV is a multi-criteria decision-making technique that ranks alternatives using concordance and discordance principles with progressive elimination.
}
\usage{
electre4 (dmatrix, bcvec, weights, p = NULL, q = NULL, v = NULL, 
    thr = 0.5, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A decision matrix where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{A vector specifying the benefit (\code{1}) or cost (\code{-1}) nature of each criterion.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{p}{Pseudo-criterion preference threshold values. Default is calculated dynamically.}
  \item{q}{Pseudo-criterion indifference threshold values. Default is set to half of \code{pt}.}
  \item{v}{Pseudo-criterion veto threshold values. Default is set to \code{pt} multiplied by 1.5.}
  \item{thr}{Threshold for credibility-based elimination. A higher threshold (e.g., 0.7 or 0.8) will make selection more restrictive, favoring only alternatives with strong credibility. A lower threshold (e.g., 0.4 or 0.5) will allow more alternatives to remain in consideration. Default threshold is \code{0.5}. }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
The ELECTRE IV method constructs concordance and discordance matrices based on pseudo-criteria thresholds. A credibility matrix is then derived, followed by a threshold-based filtering process to eliminate weaker alternatives.
}
\value{
A list containing:
\item{Concordance}{Concordance matrix representing the degree of dominance of alternatives.}
\item{Discordance}{Discordance matrix indicating veto situations.}
\item{Credibility}{Credibility matrix integrating concordance and discordance effects.}
\item{Net_Flow}{Computed net credibility flow values for ranking.}
\item{rank}{Final ranking of alternatives after threshold-based elimination.}
}
\author{
  Cagatay Cebeci
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\examples{
dmat <- matrix(c(5, 150, 3, 200, 4, 180), nrow = 3, byrow = TRUE)
bc <- c(-1, 1)
uw <- c(0.6, 0.4)

# Pseudo-criterion threshold values
pt <- list(list(p = 2), list(p = 30))
qt <- list(list(q = 1), list(q = 15))
vt <- list(list(v = 3), list(v = 45))
thr <- 0.6

reselectre4_1 <- electre4(dmatrix = dmat, bcvec = bc, weights = uw,
   p = pt, q = qt, v = vt, thr=thr)
print(reselectre4_1)

criwei <- calcweights(dmat, bcvec = bc, type = "critic")
reselectre4_2 <- electre4(dmatrix = dmat, bcvec = bc, weights = criwei)
print(reselectre4_2)
}
\concept{decision making}

\keyword{ELECTRE IV}
\keyword{decision making}
