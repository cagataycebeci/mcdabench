\name{mavt}
\alias{mavt}
\encoding{UTF-8}
\title{MAVT: Multi-Attribute Value Theory}
\description{
  Implements the Multi-Attribute Value Theory (MAVT) method for multi-criteria decision analysis. This method evaluates alternatives based on weighted aggregation of their values across multiple criteria.
}
\usage{
mavt(dmatrix, bcvec, weights, valfuncs = NULL, normvals = TRUE, 
    ss = 1, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}.}
  \item{weights}{A numeric vector of weights for each criterion. The length must equal the number of columns in \code{dmatrix}.}
  \item{valfuncs}{An optional character vector specifying the utility function for each criterion. If \code{NULL} (default), linear utility functions are used. If a single character string is provided, 
   it is applied to all criteria. Currently, \code{"linear"}, \code{"exp"}, \code{"log"}, \code{"ln"}, \code{"quad"}, and \code{"step"} are supported.}
  \item{ss}{Step size used by step function. Default is 1.}
  \item{normvals}{A logical value indicating whether to normalize the value matrix to a 0-1 range after applying the initial value functions (default is \code{TRUE}).}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
  The MAVT method involves defining a value function for each criterion that reflects the decision maker's preferences or judgments about the desirability of different levels of performance. These value functions transform the raw performance values into value scores. The overall value of each alternative is then calculated as the weighted sum of its value scores across all criteria. Alternatives are ranked based on their overall value scores in descending order.
}
\value{
  A list containing the following components:
  \item{Value_Matrix}{A numeric matrix containing the value scores of each alternative for each criterion.}
  \item{Weighted_Sum_Scores}{A numeric vector of the overall value scores for each alternative.}
  \item{Ranking}{A data frame with the ranked alternatives, their MAVT scores, and their ranks.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower rank is better).}
}
\references{
Keeney, R. L., & Raiffa, H. (1993). Decisions with multiple objectives: preferences and value trade-offs. Cambridge University Press.
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

# Apply MAVT with default linear value functions
resmavt_linear <- mavt(dmatrix=dmat, bcvec=bc, weights=uw)
print(resmavt_linear)

# Apply MAVT without normalizing the value matrix
resmavt_nonorm <- mavt(dmatrix=dmat, bcvec=bc, weights=uw, normvals = FALSE)
print(resmavt_nonorm)

}