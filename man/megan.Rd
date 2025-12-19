\name{megan}
\alias{megan}
\encoding{UTF-8}
\title{MEGAN: Multi-criteria Evaluation via Gini-weighted Attributes in Numerical Decision Matrices}
\description{
  A novel method to make multi-criteria decision analysis (MCDA) using various criteria, weights, and thresholds on numerical decision matrices.
}
\usage{
megan(dmatrix, bcvec, weights="gini", normethod = NULL, thr = NULL, 
      tht="p25", tiesmethod="average")
}
\arguments{
  \item{dmatrix}{A numerical decision matrix (\code{n x m}), where \code{n} represents alternatives and \code{m} represents criteria.}
  \item{bcvec}{A boundary condition vector, where \code{-1} represents cost criteria and \code{1} represents benefit criteria.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"gini"} because the algorithm original applies GINI weighting if no user-defined weights supplied. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be used in calculation of weights. Default is \code{NULL} because this method natively applies MaxMin normalization on distance matrix. See \code{\link{calcnormal}} for details. }
  \item{thr}{A numeric value between 0 and 100, representing the threshold for ranking alternatives based on their superiority scores. If \code{NULL}, the threshold is determined by \code{tht}. If threshold value is zero (\code{thr=0},
   the rank of alternatives is returned using \code{rank} function with the ties method \code{"average"}.}
  \item{tht}{A character string specifying the method for determining the threshold when \code{thr} is \code{NULL}. The options are \code{"sdev"} (standard deviation), \code{"p5"} (5th percentile), and \code{"p25"} (25th percentile) 
   of differences between the superiority scores of alternatives. Default is \code{"p25"}.}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The MEGAN method involves the following steps:
\enumerate{
    \item Normalize the decision matrix.
    \item Calculate criteria weights (if not provided).
    \item Calculate the weighted normalized decision matrix.
    \item Convert cost criteria to benefit criteria.
    \item Calculate distances to the minimum value for each criterion.
    \item Normalize the distance matrix.
    \item Calculate row sums of the normalized distances.
    \item Calculate superiority scores.
    \item Rank the alternatives based on their superiority scores and the given threshold. If threshold value is zero (\code{thr=0}, the rank of scores is returned using \code{rank} function with the ties method \code{"average"}.)
}
}
\value{
A list containing:
\item{dmatrix}{The original decision matrix.}
\item{wmatrix}{The weighted normalized decision matrix.}
\item{minvec}{The vector of minimum values of criteria in generalized normalized decision matrix.}
\item{distmatrix}{The matrix of distances to the minimum values.}
\item{ndistmatrix}{The normalized distance matrix.}
\item{distrowsums}{The vector of row sums of the normalized distances.}
\item{superiority}{The vector of superiority scores (as percentage over minimum) for each alternative.}
\item{rank}{The vector of ranks for each alternative.}
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\author{
  Cagatay Cebeci
}
\examples{
data(egrids)
dmat <- as.matrix(egrids$dmat) # Decision matrix
bc <- egrids$bcvec # Benefit and cost vector
wei <- egrids$weights # Criteria weights

# Apply MEGAN without threshold
resmegan_1 <- megan(dmatrix=dmat, bcvec=bc, weights=wei, thr=0)
print(resmegan_1)

# Apply MEGAN with a fixed threshold of 10%
resmegan_2 <- megan(dmatrix=dmat, bcvec=bc, weights=wei, thr=10)
print(resmegan_2$superiority)
print(resmegan_2$rank)

# Apply MEGAN with the 'p5' threshold
resmegan_3 <- megan(dmatrix=dmat, bcvec=bc, weights=wei, thr=NULL, tht="p5")
print(resmegan_3$superiority)
print(resmegan_3$rank)

# Apply MEGAN with the standard deviation threshold
resmegan_4 <- megan(dmatrix=dmat, bcvec=bc, weights=wei, thr=NULL, tht="sdev")
print(resmegan_4$superiority)
print(resmegan_4$rank)
}
