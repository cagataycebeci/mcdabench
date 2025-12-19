\name{promethee3}
\alias{promethee3}
\encoding{UTF-8}
\title{PROMETHEE III: Preference Ranking Organization METHod for Enrichment Evaluations (Method III)}
\description{
  Performs the PROMETHEE III method for complete ranking of alternatives in a Multi-Criteria Decision Analysis (MCDM) problem, allowing for interval-based ranking.
}
\usage{
promethee3(dmatrix, bcvec, weights, normethod = NULL, 
   prefuncs = NULL, thr = NULL, strict = FALSE, tiesmethod="average")
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
    An optional character string specifying the normalization method for the decision matrix because there is usually no need to normalize the decision matrix with PROMETHEE III. If not \code{NULL}, the \code{calcnormal} function will be used. Common options might include "linear", "vector", and "maxmin". Defaults to \code{NULL}, no normalization is performed. See \code{\link{calcnormal}} for details.
  }
  \item{prefuncs}{
    A character vector specifying the preference function to be used for each criterion. The length of \code{prefuncs} must be equal to the number of columns in \code{dmatrix}, or it can be a single character string to apply the same preference function to all criteria. Available preference functions are: "linear", "usual", "quasi-criterion", "v-shape", and "level".
  }
  \item{thr}{
    A list of threshold values for the preference functions. The structure depends on the chosen preference function for each criterion. It should be a list with the same length as the number of criteria. Each element of the list can be a named list or a numeric vector containing the required threshold parameters (e.g., \code{list(p = 50)} for "linear" and "v-shape", \code{list(s = 10)} for "quasi-criterion", \code{list(q = 5, p = 20)} for "level"). If a single list or numeric vector is provided, it will be applied to all criteria. If \code{NULL}, default thresholds based on the maximum absolute difference for each criterion will be used (with a warning).
  }
  \item{strict}{
     A logical value indicating whether to perform strict ranking (\code{TRUE}) or interval-based ranking (\code{FALSE}). If \code{TRUE}, the result is identical to PROMETHEE II. Default is \code{FALSE}.
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
  The PROMETHEE III method extends PROMETHEE II to allow for interval-based ranking of alternatives, providing a more nuanced result that can account for uncertainty.  It calculates positive (Phi_plus), negative (Phi_minus), and net (Phi_net) outranking flows. The final ranking depends on the \code{strict} parameter.  If \code{strict = FALSE}, alternatives are ranked based on intervals (which, in this basic implementation, are point estimates equal to Phi_net). If \code{strict = TRUE}, the ranking is based on Phi_net, as in PROMETHEE II.
}
\value{
  A list containing the following components:
  \item{Phi_plus}{A numeric vector of positive outranking flows for each alternative.}
  \item{Phi_minus}{A numeric vector of negative outranking flows for each alternative.}
  \item{Phi_net}{A numeric vector of net outranking flows for each alternative.}
  \item{rank}{
    A numeric vector representing the complete ranking of alternatives. If \code{strict = FALSE}, this represents an interval-based ranking (in this implementation, the intervals are point estimates). If \code{strict = TRUE}, this is a strict ranking based on \code{Phi_net} values (identical to PROMETHEE II). Ties are handled using the "average" method.
    }
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
  \code{\link{promethee1}}, \code{\link{promethee2}}, \code{\link{promethee4}}, \code{\link{promethee5}}, \code{\link{promethee6}}
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

# User-defined weights
uw <- c(0.6, 0.4)

# Threshold values
tvals <- list(list(p = 50), list(p = 30))

# Preference functions
pf <- c("linear", "linear")

# Run PROMETHEE III function (strict = FALSE for interval ranking)
respromethee3_interval <- promethee3(dmatrix=dmat, bcvec=bc, weights=uw, 
   prefuncs = pf, thr = tvals, strict = FALSE)
print(respromethee3_interval)

# Run PROMETHEE III function (strict = TRUE for strict ranking)
respromethee3_strict <- promethee3(dmatrix=dmat, bcvec=bc, weights=uw, 
   prefuncs = pf, thr = tvals, strict = TRUE)
print(respromethee3_strict)
}
