\name{promethee1}
\alias{promethee1}
\encoding{UTF-8}
\title{PROMETHEE I Method for Partial Ranking in Multi-Criteria Decision Analysis}
\description{
Implements the PROMETHEE I method to compute positive and negative preference flows, generate pairwise outranking comparisons, and determine a partial ranking of alternatives.
}
\usage{
promethee1(dmatrix, bcvec, weights, normethod = NULL, 
    prefuncs = NULL, thr = NULL, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric matrix or data frame representing the decision matrix. Rows are alternatives; columns are criteria.}
  \item{bcvec}{A numeric vector of \code{-1} (cost) or \code{1} (benefit) values specifying the optimization direction of each criterion.}
  \item{weights}{A numeric vector of weights for criteria. Must sum to approximately 1.}
  \item{normethod}{Optional normalization method passed to \code{\link{calcnormal}}. If \code{NULL}, no normalization is applied.}
  \item{prefuncs}{A character vector or single value specifying preference function types for each criterion. Supported types include
     \code{"linear"}, \code{"usual"}, \code{"quasi"}, \code{"v-shape"}, \code{"level"}, \code{"linear-indif"}, and \code{"gaussian"}.}
  \item{thr}{A list of threshold parameters for preference functions. Each element may include \code{p}, \code{q}, \code{s}, or \code{sigma} 
     depending on the function type.}
  \item{tiesmethod}{Character string specifying tie-breaking method for final Phi+ ranking. Passed to \code{\link[base]{rank}}.}
}
\value{
A list containing:
  \item{Phi_plus}{A vector of positive preference flow values for each alternative.}
  \item{Phi_minus}{A vector of negative preference flow values for each alternative.}
  \item{Partial_Ranking}{A matrix indicating pairwise outranking relations (\code{"P"}, \code{"Q"}, \code{"I"}, \code{"R"}, \code{"="}).}
  \item{rank}{A numeric vector of ranks based on descending Phi+ values.}
}
\details{
PROMETHEE I allows pairwise comparison of alternatives under multiple criteria using weighted preference functions. The method computes how much an alternative is preferred to others (Phi+) and how much it is outranked (Phi-), without producing a complete ordering in all cases. Indifference and incomparability are explicitly represented.

You can customize the behavior via different preference functions and thresholds for each criterion. When no thresholds are provided, default ones based on column differences are used.
}
\examples{
# Sample decision matrix
dm <- matrix(c(10, 5, 8,
               6, 9, 7,
               8, 7, 9), nrow = 3, byrow = TRUE)
rownames(dm) <- paste0("A", 1:3)
colnames(dm) <- paste0("C", 1:3)

# All criteria are benefits
bc <- c(1, 1, 1)

# Weights sum to 1
userwei <- c(0.3, 0.4, 0.3)

# Preference functions and thresholds
prefuncs <- c("linear", "v-shape", "level")
thr <- list(
  list(p = 5),         # Linear
  list(p = 4),         # V-shape
  list(q = 1, p = 3)   # Level
)

# Apply PROMETHEE I
respromethee1 <- promethee1(dmatrix=dm, bcvec=bc, weights=userwei, 
   prefuncs = prefuncs, thr = thr)
print(respromethee1)
}
\author{
Cagatay Cebeci
}
