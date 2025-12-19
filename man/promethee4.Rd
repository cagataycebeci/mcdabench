\name{promethee4}
\alias{promethee4}
\encoding{UTF-8}
\title{PROMETHEE IV: Preference Ranking Organization METHod for Enrichment Evaluations (Method IV)}
\description{
  Performs the PROMETHEE IV method for complete ranking of alternatives in a Multi-Criteria Decision Analysis (MCDM) problem, considering the average outranking flow.
}
\usage{
promethee4(dmatrix, bcvec, weights, normethod = NULL, 
   prefuncs = NULL, thr = NULL, alpha = 0.2, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{
    A numeric matrix representing the decision matrix, where rows correspond to alternatives and columns correspond to criteria.
  }
  \item{bcvec}{
    A numeric vector indicating the benefit or cost criterion. Use 1 for benefit (higher is better) and -1 for cost (lower is better) for each criterion. The length of \code{bcvec} must be equal to the number of columns in \code{dmatrix}.
  }
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{normethod}{
    An optional character string specifying the normalization method for the decision matrix because there is usually no need to normalize the decision matrix with PROMETHEE IV. If not \code{NULL}, the \code{calcnormal} function will be used. Common options might include "linear", "vector", etc. Defaults to \code{NULL}, no normalization is performed. See \code{\link{calcnormal}} for details.
  }
  \item{prefuncs}{
    A character vector specifying the preference function to be used for each criterion. The length of \code{pref_functions} must be equal to the number of columns in \code{dmatrix}, or it can be a single character string to apply the same preference function to all criteria. Available preference functions are: "linear", "usual", "quasi-criterion", "v-shape", and "level".
  }
  \item{thr}{
    A list of threshold values for the preference functions. The structure depends on the chosen preference function for each criterion. It should be a list with the same length as the number of criteria. Each element of the list can be a named list or a numeric vector containing the required threshold parameters (e.g., \code{list(p = 50)} for "linear" and "v-shape", \code{list(s = 10)} for "quasi-criterion", \code{list(q = 5, p = 20)} for "level"). If a single list or numeric vector is provided, it will be applied to all criteria. If \code{NULL}, default thresholds based on the maximum absolute difference for each criterion will be used (with a warning).
  }
  \item{alpha}{
    A numeric value between 0 and 1 representing the influence of the average outranking flow in the final ranking. A higher value gives more importance to the average flow. Default value is 0.2.
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
  PROMETHEE IV ranks alternatives by considering their deviation from the average outranking flow.  It calculates positive (\code{Phi_plus}), negative (\code{Phi_minus}), and net (\code{Phi_net}) outranking flows, and then computes the average net outranking flow (\code{Phi_avg}). The final ranking is based on the adjusted net outranking flow: \code{Phi_net - alpha * Phi_avg}.
}
\value{
  A list containing the following components:
  \item{Phi_plus}{A numeric vector of positive outranking flows for each alternative.}
  \item{Phi_minus}{A numeric vector of negative outranking flows for each alternative.}
  \item{Phi_net}{A numeric vector of net outranking flows for each alternative.}
  \item{Phi_avg}{The average of the net outranking flows.}
  \item{rank}{A numeric vector representing the complete ranking of alternatives, considering the average outranking flow.}
}
\references{
Brans, J. P., & Mareschal, B. (1988). Geometrical representations for MCDA. \emph{European Journal of Operational Research, 34}(1), 69-77. <doi:10.1016/0377-2217(88)90456-0>

Brans, J. P., Mareschal, B. (2005). Promethee Methods. In: \emph{Multiple Criteria Decision Analysis: State of the Art Surveys. 
 International Series in Operations Research & Management Science, vol 78.} Springer, New York, NY. <doi:10.1007/0-387-23081-5_5>

Behzadian, M., Kazemzadeh, R. B., Albadvi, A., & Aghdasi, M. (2010). PROMETHEE: A comprehensive literature review on methodologies and applications.
 \emph{European Journal of Operational Research, 200}(1), 198-215. <doi:10.1016/j.ejor.2009.01.021>

Brans, J. P., De Smet, Y. (2016). PROMETHEE Methods. In: \emph{Greco, S., Ehrgott, M., Figueira, J. (eds) Multiple Criteria Decision Analysis. 
 International Series in Operations Research & Management Science, 233}. Springer, New York, NY. <doi:10.1007/978-1-4939-3094-4_6>
}
\seealso{
  \code{\link{promethee1}}, \code{\link{promethee2}}, \code{\link{promethee3}}, \code{\link{promethee5}}, \code{\link{promethee6}}
}
\author{
 Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(5, 150, 3, 200, 4, 180), nrow = 3, byrow = TRUE)
rownames(dmat) <- paste0("A", 1:3)
colnames(dmat) <- c("C1", "C2")

# Benefit-Cost vector
bc <- c(-1, 1)

# Weights
uw <- c(0.6, 0.4)

# Threshold values
tvals <- list(list(p = 50), list(p = 30))

# Preference function
pf <- c("linear", "linear")

# Run PROMETHEE IV function with alpha = 0.2
respromethee4_1 <- promethee4(dmatrix=dmat, bcvec=bc, weights=uw,
  thr= tvals, prefuncs = pf, alpha = 0.2)
print(respromethee4_1)

# Run PROMETHEE IV function with alpha = 0.8
respromethee4_2 <- promethee4(dmatrix=dmat, bcvec=bc, weights=uw,
  thr = tvals, prefuncs = pf, alpha = 0.8)
print(respromethee4_2)
}
