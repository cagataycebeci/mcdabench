\name{waspas}
\alias{waspas}
\encoding{UTF-8}
\title{WASPAS: Weighted Aggregated Sum Product Assessment}
\description{
An implementation of the Weighted Aggregated Sum Product Assessment (WASPAS) method, which combines the Weighted Sum Model (WSM) and the Weighted Product Model (WPM), to evaluate and rank alternatives based on normalized scores and assigned weights.
}
\usage{
waspas(dmatrix, bcvec, weights, normethod = "linear", v = 0.5, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{
    Decision matrix where rows represent alternatives and columns represent criteria.
  }
  \item{bcvec}{
    Benefit-cost vector, specifying for each criterion whether it is a benefit (\code{1}) or a cost (\code{-1}).
  }
  \item{weights}{A numeric vector of user-defined weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method to be applied. Defaults to \code{"linear"}. See \code{\link{calcnormal}} for details. }
  \item{v}{
    A parameter that balances the contributions of WSM and WPM. A value between 0 and 1, where \code{v = 0.5} gives equal weight to both models. Default is \code{0.5}.
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The WASPAS method integrates the advantages of the Weighted Sum Model (WSM) and the Weighted Product Model (WPM) through a parameterized combination. Alternatives are evaluated based on a sum-based score (\code{qsum}) and a product-based score (\code{qprod}), with the final score (\code{p}) calculated as a weighted average of the two.
}
\value{
A list containing the following components:
  \item{nmatrix}{The normalized decision matrix.}
  \item{qsum}{The aggregated score based on the Weighted Sum Model.}
  \item{qprod}{The aggregated score based on the Weighted Product Model.}
  \item{p}{The final score for each alternative, calculated as a combination of \code{qsum} and \code{qprod}.}
 \item{rank}{The ranks of final scores for alternatives.}
}
\references{
Chakraborty, S., Zavadskas, E. K., & Antuchevičienė, J. (2015). Applications of WASPAS method as a multi-criteria decision-making tool. URI <https://etalpykla.vilniustech.lt/handle/123456789/151097>.

Zavadskas, E. K., Vilutienė, T., Turskis, Z., & Šaparauskas, J. (2014). Multi-criteria analysis of Projects’ performance in construction. 
 \emph{Archives of Civil and Mechanical Engineering, 14}, 114-121. <doi:10.1016/j.acme.2013.07.006>
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

# Apply WASPAS with user-defined weights
reswaspas_1 <- waspas(dmatrix = dmat, bcvec=bc, 
   weights=uw, normethod = "linear", v = 0.5)
print(reswaspas_1)

# Apply WASPAS with GINI weights
giniw <- calcweights(dmat, bcvec=bc, type="gini")
reswaspas_2 <- waspas(dmatrix = dmat, bcvec=bc, 
   weights=giniw, normethod = "linear", v = 0.5)
print(reswaspas_2)
}
