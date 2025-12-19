\name{maut}
\alias{maut}
\encoding{UTF-8}
\title{MAUT: Multi-Attribute Utility Theory}
\description{
  The MAUT evaluates alternatives based on weighted aggregation of their utilities across multiple criteria.
}
\usage{
maut(dmatrix, bcvec, weights, utilfuncs = NULL, normutil = TRUE,  
    ss = 1, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric matrix representing the decision matrix. Rows are alternatives and columns are criteria.}
  \item{bcvec}{A numeric vector indicating the benefit/cost nature of each criterion. Use 1 for benefit (to be maximized) and -1 for cost (to be minimized). The length of \code{bcvec} must equal the number of columns in \code{dmatrix}.}
  \item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"}. See \code{\link{calcweights}} for details.}
  \item{utilfuncs}{An optional character vector specifying the utility function for each criterion. If \code{NULL} (default), linear utility functions are used. If a single character string is provided, 
   it is applied to all criteria. Currently, \code{"linear"}, \code{"exp"}, \code{"log"}, \code{"ln"}, \code{"quad"}, and \code{"step"} are supported.}
  \item{ss}{Step size used by step utility function. Default is 1.}
  \item{normutil}{A logical value indicating whether to normalize the utility matrix to a 0-1 range after applying the initial utility functions (default is \code{TRUE}).}
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.}
}
\details{
 The MAUT method involves defining a utility function for each criterion that reflects the decision maker's preferences. These utility functions transform the raw performance values into utility scores. The overall 
 utility of each alternative is then calculated as the weighted sum of its utility scores across all criteria. Alternatives are ranked based on their overall utility scores in descending order.

 In MAUT based on utility theory, the decision-maker calculates the utility of each criterion and strives to best reflect their preferences in the decision-making process. Utility is shaped by subjective perceptions 
 and risk factors. On the other hand, MAVT based on value theory, the decision-maker determines the value of each criterion and optimizes these values to achieve specific objectives. MAVT generally focuses on the direct
 value of each criterion and adopts a more objective approach. Conceptual differences are related to the decision-maker’s preferences and risk perceptions. While working with the same matrix, different decision-makers may 
 use either MAUT or MAVT, but the chosen approach affects how criteria should be evaluated. For example, if a decision-maker is risk-averse, they may prefer to make decisions using MAUT, as it flexibly models utility functions 
 by taking into account risk aversion or risk-taking. On the other hand, if a decision-maker prefers a more objective approach and wants to directly optimize the value of each criterion, MAVT may be more suitable. In this case, 
 the obtained value for each criterion takes precedence.
}
\value{
  A list containing the following components:
  \item{Utility_Matrix}{A numeric matrix containing the utility scores of each alternative for each criterion.}
  \item{Weighted_Sum_Scores}{A numeric vector of the overall utility scores for each alternative.}
  \item{rank}{A numeric vector representing the rank of each alternative (lower rank is better).}
}
\references{
Keeney, R. L., & Raiffa, H. (1993). Decisions with multiple objectives: preferences and value trade-offs. Cambridge University Press.
Ishizaka, A., & Nemery, P. (2013). Multi-criteria decision analysis: methods and software. John Wiley & Sons.
}
\seealso{
  \code{\link{calcweights}}
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

# Utility functions
uf <- c("linear", "linear","quad","step")

# Apply MAUT with default linear utility functions
resmaut_linear <- maut(dmatrix=dmat, bcvec=bc, weights=uw)
print(resmaut_linear)

# Apply MAUT without normalizing the utility matrix
resmaut_nonorm <- maut(dmatrix=dmat, bcvec=bc, weights=uw, normutil = FALSE)
print(resmaut_nonorm)

# Apply MAUT without normalizing the utility matrix with utility functions
resmaut_utils <- maut(dmatrix=dmat, bcvec=bc, weights=uw, utilfuncs=uf, normutil = FALSE)
print(resmaut_utils)
}