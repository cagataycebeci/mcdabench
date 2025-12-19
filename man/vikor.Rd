\name{vikor}
\alias{vikor}
\encoding{UTF-8}
\title{VIKOR: VIsekriterijumska optimizacija i KOmpromisno Resenje}
\description{
Implements the VIsekriterijumska optimizacija i KOmpromisno Resenje (VIKOR) (Multi-criteria optimization and Compromise Solution in English) method to rank alternatives by evaluating their compromise solutions based on the closeness to ideal solutions.
}
\usage{
vikor(dmatrix, bcvec, weights, normethod = "maxmin", v = 0.5, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{
    Decision matrix where rows represent alternatives and columns represent criteria.
  }
  \item{bcvec}{
    Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).
  }
  \item{weights}{A numeric vector of user-defined weights representing the relative importance of each criterion, or a character string to calculate it internally. Defaults to \code{"equal"} to apply equal weighting for criteria. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"maxmin"}. See \code{\link{calcnormal}} for details. }
  \item{v}{
    The weight of the strategy for maximum group utility. A value between 0 and 1. Default is \code{0.5}.
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The VIKOR method normalizes the decision matrix, applies weights, and calculates the utility measure (\code{S}), the regret measure (\code{R}), and the compromise measure (\code{Q}). These measures are used to rank the alternatives, with the best compromise solution identified as the one with the smallest \code{Q} value.
}
\value{
A list containing the following components:
  \item{\code{nmatrix}}{The normalized decision matrix.}
  \item{\code{fminus}}{A vector of the worst performance values (minimum values) for each criterion.}
  \item{\code{fstar}}{A vector of the best performance values (maximum values) for each criterion.}
  \item{\code{S}}{A vector of the group utility values for each alternative.}
  \item{\code{R}}{A vector of the individual regret values for each alternative.}
  \item{\code{Q}}{A vector of the compromise values for each alternative.}
  \item{\code{rank}}{The ranking of alternatives based on the Q values.}
}
\references{
Opricović, S. (1979). Multicriteria Optimization in Civil Engineering. (Ph.D. Dissertation, in Serbian), Faculty of Civil Engineering, Belgrade.

Duckstein, L., & Opricovic, S. (1980). Multiobjective optimization in river basin development. \emph{Water Resources Research, 16}(1), 14-20.

Opricović, S. (1998). Multicriteria optimization of civil engineering systems. \emph{Faculty of civil engineering, Belgrade, 2}(1), 5-21.

Opricovic, S., & Tzeng, G. H. (2004). Compromise solution by MCDM methods: A comparative analysis of VIKOR and TOPSIS. 
 \emph{European journal of Operational Research, 156}(2), 445-455. <doi:10.1016/S0377-2217(03)00020-1>
}
\seealso{
  \code{\link{calcweights}}, \code{\link{calcnormal}}
}
\author{
  Cagatay Cebeci
}
\examples{
# Sample decision matrix
dmat <- matrix(c(
  84, 10, 5, 20,  
  66, 18, 8, 15,
  74, 12, 6, 25,
  90, 22, 9, 18, 
  68, 18, 7, 15  
), nrow = 5, byrow = TRUE)
rownames(dmat) <- c("A1", "A2", "A3", "A4", "A5")
colnames(dmat) <- c("C1", "C2", "C3", "C4")

# Benefit-cost vector (1: benefit, -1: cost)
bc <- c(1, -1, 1, -1)

# User-defined weights
uw <- c(0.1, 0.4, 0.3, 0.2)

resvikor_1 <- vikor(dmatrix=dmat, bcvec=bc, weights=uw, v = 0.5)
print(resvikor_1)

resvikor_2 <- vikor(dmatrix=dmat, bcvec=bc, weights=uw, v = 0.1)
print(resvikor_2)

criw <- calcweights(dmat, bcvec=bc, type="critic")
resvikor_3 <- vikor(dmatrix=dmat, bcvec=bc, weights=criw, 
  normethod="vector", v = 0.5)
print(resvikor_3)
}
