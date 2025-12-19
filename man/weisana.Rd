\name{weisana}
\alias{weisana}
\encoding{UTF-8}
\title{Weight Sensitivity Analysis for Multi-Criteria Decision-Making}
\description{
Performs a weight sensitivity analysis by generating modified weight scenarios
and applying a user-defined MCDA ranking method (e.g., TOPSIS, VIKOR) to each.
Summarizes unique ranking patterns in rankings.
}
\usage{
weisana(dmatrix, bcvec, weights, weimethod = "gradual", weipars = NULL,
        mcdamethod = topsis, methodpars, sensplot = TRUE)
}
\arguments{
  \item{dmatrix}{Decision matrix (numeric matrix or data frame) where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{Numeric vector with values -1 or 1, indicating whether each criterion is a cost (-1) or benefit (+1).}
  \item{weights}{Numeric vector of weights for criteria. If omitted, equal weights are assigned.}
  \item{weimethod}{Weight perturbation method: \code{"gradual"} or \code{"random"}.}
  \item{weipars}{List of parameters for the weight generation method. Includes \code{rp}, \code{ss}, and \code{niters}.}
  \item{mcdamethod}{Ranking function to apply, such as \code{topsis} or \code{vikor}. Must return a list with a \code{$rank} vector.}
  \item{methodpars}{Named list of additional parameters passed to \code{mcdamethod}.}
  \item{sensplot}{Plot for ranking distributions obtained with different weighting scenarios.}
}
\value{
  A named list with the following components:
  \item{computing_time}{Time taken to perform ranking computations.}
  \item{weights_matrix}{Matrix of generated weight scenarios.}
  \item{ranking_matrix}{Matrix of rankings per weight set.}
  \item{sensitivity_table}{Data frame showing unique ranking patterns, counts, and percentages.}
}
\details{
This function allows for analyzing the robustness of MCDA results under different weight variations.
It supports both gradual and random perturbation schemes and is agnostic to the MCDA method used.
}
\examples{
# Sample decision matrix and parameters
dmat <- matrix(runif(20), nrow=5)
rownames(dmat) <- paste0("A", 1:nrow(dmat))
colnames(dmat) <- paste0("C", 1:ncol(dmat))
print(dmat)
# Benefit-cost vector
bc <- c(1, 1, -1, 1)

# A very simple custom MCDA function
customfunc <- function(dmatrix, bcvec, weights, v=0.5, tiesmethod="average") {
  scores <- (rowSums(scale(dmatrix) * weights))^v
  ranks <- rank(-scores, ties.method=tiesmethod)
  list(rank = ranks)
}

# Sensitivity of the custom MCDM with gradual weight modification
mp <- list(v = 0.6, tiesmethod="average")
wp <- list(rp = c(0.1, 0.5, 0.7))
resgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = customfunc, methodpars = mp)
print(resgrawei)

# Sensitivity of the custom MCDM with random weight modification
mp <- list(v = 0.6)
wp <- list(ss = 0.1, niters=20)
resrandwei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "random", weipars = wp,
     mcdamethod = customfunc, methodpars = mp)
print(resrandwei)

# Test TOPSIS
mp <- list()
wp <- list(rp = seq(0.01, 0.5, 0.05))

topsisgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = topsis, methodpars = mp)
print(topsisgrawei)

# Test VIKOR
mp <- list(v=0.5)
wp <- list(rp = seq(0.01, 0.5, 0.05))

vikorgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = vikor, methodpars = mp, sensplot=FALSE)
print(vikorgrawei)
sensplot(vikorgrawei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis for VIKOR", colpal=terrain.colors(10))

# Test WASPAS
mp <- list(v=0.5, normethod="linear", tiesmethod="average")
wp <- list(rp = seq(0.01, 0.6, 0.01))

waspasgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = waspas, methodpars = mp)
print(waspasgrawei)
sensplot(waspasgrawei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis for WASPAS", colpal=terrain.colors(10))
   
waspasrandwei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "random", weipars = list(ss=0.05, niters=50),
     mcdamethod = waspas, methodpars = mp)
print(waspasrandwei)
sensplot(waspasrandwei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis (random) for WASPAS", colpal=terrain.colors(10))
}

\author{
Cagatay Cebeci
}
