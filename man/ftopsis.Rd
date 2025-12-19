\name{ftopsis}
\alias{ftopsis}
\encoding{UTF-8}
\title{Fuzzy TOPSIS Method}
\description{
Computes rankings for decision alternatives using the Fuzzy TOPSIS method, which evaluates alternatives based on triangular fuzzy numbers (l, m, u).
}
\usage{
ftopsis(fuzzydmatrix, bcvec, fuzzyweights, tiesmethod = "average")
}
\arguments{
  \item{fuzzydmatrix}{A decision matrix where each criterion is represented by triangular fuzzy numbers (l, m, u), 
   structured as three adjacent columns per criterion.}
  \item{bcvec}{A vector indicating the type of each criterion: 1 for benefit criteria and -1 for cost criteria. 
   Must have the same length as the number of criteria.}
  \item{fuzzyweights}{A matrix where each row corresponds to a criterion and contains three columns representing fuzzy weights (l, m, u).}
  \item{tiesmethod}{Method for resolving ties in ranking, defaulting to `"average"`.}
}
\value{
A list containing:
  \item{Criteria}{Number of criteria evaluated.}
  \item{fuzzyweights}{Fuzzy weights applied to each criterion.}
  \item{FPIS}{Fuzzy Positive Ideal Solution.}
  \item{FNIS}{Fuzzy Negative Ideal Solution.}
  \item{distance_to_FPIS}{Euclidean distance of each alternative to FPIS.}
  \item{distance_to_FNIS}{Euclidean distance of each alternative to FNIS.}
  \item{closeness_coefficient_fuzzy}{Closeness coefficient values for ranking.}
  \item{rank}{Final ranking of alternatives based on closeness coefficients.}
}
\references{
Nadaban, S., Dzitac, S., & Dzitac, I. (2016). Fuzzy TOPSIS: a general view. \emph{Procedia Computer Science, 91}, 823-831. <doi:10.1016/j.procs.2016.07.088>
}
\examples{
fuzzydmatrix <- matrix(c(
  8, 10, 12,    # A1-C1
  18, 20, 22,   # A2-C1
  28, 30, 32,   # A3-C1
  38, 40, 42,   # A4-C1
  45, 50, 55,   # A1-C2
  55, 60, 65,   # A2-C2
  35, 40, 45,   # A3-C2
  25, 30, 35,   # A4-C2
  28, 30, 32,   # A1-C3
  18, 20, 22,   # A2-C3
  38, 40, 42,   # A3-C3
  48, 50, 52,   # A4-C3
  75, 80, 85,   # A1-C4
  65, 70, 75,   # A2-C4
  55, 60, 65,   # A3-C4
  45, 50, 55    # A4-C4
), nrow = 4, byrow = TRUE)

bcvec <- c(1, 1, -1, 1) 

fuzzyweights <- matrix(c(
  0.25, 0.30, 0.35, 
  0.35, 0.40, 0.45, 
  0.05, 0.10, 0.15, 
  0.15, 0.20, 0.25  
), nrow = 4, byrow = TRUE)

resftopsis <- ftopsis(fuzzydmatrix, bcvec, fuzzyweights)
print(resftopsis)
}

\author{
Cagatay Cebeci
}
