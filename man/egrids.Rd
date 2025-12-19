\name{egrids}
\alias{egrids}
\title{System Configuration for Optimizing Smart Grids}
\description{
The \code{egrids} dataset contains a decision matrix (\code{egrids}), a benefit-cost vector (\code{bcvec}), and a weights vector (\code{weights}) designed for analyzing and ranking grid alternatives (Gs) using multi-criteria decision-making methods.
}
\usage{
data(egrids)
}
\format{
A list with the following components:
\describe{
  \item{\code{dmat}}{
    A numeric matrix where rows represent grid alternatives (G) and columns represent criteria. Criteria include:
    \itemize{
      \item \code{C1}: Energy Efficiency (%). Higher is better
      \item \code{C2}: System Reliability (%). Higher is better
      \item \code{C3}: Renewable Integration (%). Higher is better
      \item \code{C4}: Grid Stability Score. Higher is better
      \item \code{C5}: Environmental Impact (Score). Lower is better
      \item \code{C6}: Operational Cost ($M). Lower is better
      \item \code{C7}: Emissions (\eqn{CO_2} tons/year). Lower is better
      \item \code{C8}: Response Time (seconds). Lower is better
      \item \code{C9}: Infrastructure Investment ($M). Lower is better
      \item \code{C10}: Maintenance Cost ($M). Lower is better
    }
  }
  \item{\code{bcvec}}{
    A numeric benefit and cost vector specifying the type of each criterion:
    \itemize{
      \item \code{1}: Benefit criterion meaning higher value is better.
      \item \code{-1}: Cost criterion meaning lower value is better.
    }
  }
  \item{\code{weights}}{
    A numeric vector containing the weights assigned to each criterion.
  }
}
}
\details{
This synthetic dataset provides basic criteria for evaluating and comparing smart grids using various multi-criteria decision analysis methods such as TOPSIS, VIKOR, WASPAS, MEGAN, and others. Each row in the decision matrix (\code{dmat}) corresponds to a different grid, while each column represents a specific criterion.

The \code{bcvec} vector identifies whether each criterion is a benefit (\code{1}) or a cost (\code{-1}), enabling proper handling in decision-making methods.

The \code{weights} vector contains the weights given to criteria.
}
\author{
Cagatay Cebeci
}
\examples{
data(egrids)
dmat <- egrids$dmat # Decision matrix
bcvec <- egrids$bcvec # Benefit and cost vector
weights <- egrids$weights # Criteria weights

print(dmat) 
print(bcvec) 
print(weights) 
}
\keyword{dataset}
\keyword{smart grids}
\keyword{decision-making}