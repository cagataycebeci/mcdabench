\name{mcdabench-package}
\alias{mcdabench-package}
\title{Benchmarking for Multi-Criteria Decision Analysis Methods}
\description{
Comparative benchmarking of Multi-Criteria Decision Analysis (MCDA) methods. These methods are designed to evaluate and rank alternatives based on multiple criteria, applying normalization, weighting, and aggregation techniques. The package includes popular decision-making methods such as MEGAN, TOPSIS, VIKOR, WASPAS, WSM, WPM, and others, facilitating flexible and efficient analyses for multi-criteria problems.
}
\details{
The \code{mcdabench} package aims to simplify the application of MCDA methods for researchers, analysts, and decision-makers. The package includes functionality for:

\itemize{
  \item Data normalization (e.g., \code{"maxmin"}, \code{"sum"}, \code{"linear"})
  \item Weight calculation using various techniques (\code{"equal"}, \code{"entropy"}, \code{"critic"}, etc.)
  \item Methods for ranking alternatives based on multiple criteria:
    \itemize{
      \item \code{edas} - Evaluation based on Distance from Average Solution
      \item \code{mabac} - Multi-Attributive Border Approximation Area Comparison
      \item \code{marcos} - Measurement of Alternatives and Ranking According to Compromise Solution
      \item \code{ram} - Root Assessment Method
      \item \code{topsis} - Technique for Order of Preference by Similarity to Ideal Solution
      \item \code{vikor} - Compromise ranking based on utility and regret measures
      \item \code{waspas} - Combines Weighted Sum Model (WSM) and Weighted Product Model (WPM)
      \item \code{wsm} - Weighted Sum Method
      \item \code{wpm} - Weighted Product Method
    }
  \item Flexible support for benefit-cost criteria, customizable weight schemes, and normalization methods.
}

This package is suited for applications in various decision-making contexts, such as supply chain management, project evaluation, and resource allocation.

For detailed examples and method-specific information, refer to the help pages of individual functions.
}
\author{
Cagatay Cebeci
}
\examples{
data(egrids)
dm <- as.matrix(egrids$dmat) # Decision matrix
bc <- egrids$bcvec # Benefit and cost vector
wei <- egrids$weights # Weigths of the criteria

# Perform MEGAN analysis
resmegan <- megan(dmatrix=dm, bcvec=bc, weights=wei, thr=0)
print(resmegan$superiority)  # MEGAN superiority scores
print(resmegan$rank) # Ranks of alternatives by MEGAN

# Perform VIKOR analysis
resvikor <- vikor(dmatrix=dm, bcvec=bc, weights=wei, v=0.5)
print(resvikor$Q)  # VIKOR compromise scores
print(resvikor$rank)  # Ranks of alternatives by VIKOR

# Perform TOPSIS analysis
restopsis <- topsis(dmatrix=dm, bcvec=bc, weights=wei)
print(restopsis$p)  # TOPSIS closeness coefficient
print(restopsis$rank) # Ranks of alternatives by TOPSIS
}
\concept{Multi-Criteria Decision Analysis}
\concept{Multi-Criteria Decision Making}
\keyword{MCDA}
\keyword{MCDM}
\keyword{decision-making}
\keyword{multi-criteria}