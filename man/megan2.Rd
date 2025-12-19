\name{megan2}
\alias{megan2}
\title{Multi-Criteria Evaluation with Gradually-weighted Aggregation}
\description{
  Performs multi-criteria evaluation with gradually modified weights, aggregating numerical decision matrices for ranking.
}
\usage{
megan2 (dmatrix, bcvec, weights, weimethod="random", normethod=NULL, thr=0, tht="p25", 
   tiesmethod="average", ss=0.05, niters=10, rp=c(0.1, 0.5, 0.1))
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix, where rows correspond to alternatives and columns correspond to criteria.}
  \item{bcvec}{A numeric vector indicating the benefit or cost criterion. Use 1 for benefit (higher is better) and -1 for cost (lower is better).}
  \item{weights}{A numeric vector of weights representing the relative importance of each criterion. If omitted, equal weights are assigned.}
  \item{normethod}{Normalization method applied to criteria before aggregation. If \code{NULL}, raw data is used.}
  \item{weimethod}{Weight modification method applied to chanage criteria weights. If \code{gradual}, gradual weighting method is used. Default is \code{random} for random change using a step size.}
  \item{thr}{Threshold value applied in preference aggregation. Default is 0.}
  \item{tht}{Threshold type used in preference evaluation. Defaults to "p25".}
  \item{tiesmethod}{Method for handling ties in ranking. Defaults to "average".}
  \item{ss}{Step size for modifying weights randomly. Defaults to 0.05.}
  \item{niters}{Number of iterations for generating random weight modifications. Defaults to 10.}
  \item{rp}{A numeric vector defining the percentage reductions applied to each criterion's weight for gradual modifications. Defaults to \code{c(0.1, 0.5, 0.1)}.}
}
\details{
  This function applies gradual modifications to the initial weights, redistributing the importance of criteria step by step. 
  Rankings are aggregated based on preference evaluations and adjusted to handle ties.
}
\value{
  A list containing:
  \item{weimat}{A matrix containing different weight scenarios generated through gradual modifications.}
  \item{rankmat}{A matrix of rankings corresponding to different weight scenarios.}
  \item{rank}{Final aggregated ranking of alternatives based on preference evaluations.}
}
\seealso{
  \code{\link{gengradwei}}, \code{\link{rankaggregate}}, \code{\link{megan}}
}
\author{
  Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(5, 150, 3, 200, 4, 180), nrow = 3, byrow = TRUE)
rownames(dmat) <- paste0("A", 1:3)
colnames(dmat) <- c("C1", "C2")

# Benefit-Cost vector
bc <- c(-1, 1)

# Initial weights
init_weights <- c(C1 = 0.6, C2 = 0.4)

# Run megan2 function with gradual weight modifications
result_grad <- megan2(dmatrix = dmat, bcvec = bc, weights = init_weights, rp = c(0.1, 0.3, 0.2))

# Generated weight scenarios
print(result_grad$weimat)

print(result_grad$rank) # Aggregated rankings

# Run megan2 function with random weight modifications
result_rand <- megan2(dmatrix = dmat, bcvec = bc, weights = init_weights, ss=0.02, niters=10)

# Generated weight scenarios
print(result_rand$weimat)

print(result_rand$rank) # Aggregated rankings

}

