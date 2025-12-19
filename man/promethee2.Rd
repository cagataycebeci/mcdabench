\name{promethee2}
\alias{promethee2}
\encoding{UTF-8}
\title{PROMETHEE II: Preference Ranking Organization METHod for Enrichment Evaluations (Method II)}
\description{
  Performs the PROMETHEE II method for complete ranking of alternatives in a Multi-Criteria Decision Analysis (MCDA) problem.
}
\usage{
promethee2(dmatrix, bcvec, weights, normethod = NULL, prefuncs = NULL, 
    thr = NULL, tiesmethod="average")
}
\arguments{
  \item{dmatrix}{
    A numeric matrix representing the decision matrix, where rows correspond to alternatives and columns correspond to criteria.
  }
  \item{bcvec}{
    A numeric vector indicating the benefit or cost criterion. Use 1 for benefit (higher is better) and -1 for cost (lower is better) for each criterion. The length of \code{bcvec} must be equal to the number of columns in \code{dmatrix}.
  }
  \item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
  \item{normethod}{The normalization method if applied. Defaults to \code{NULL}. See \code{\link{calcnormal}} for details. }
  \item{prefuncs}{
    A character vector specifying the preference function to be used for each criterion. The length of \code{prefuncs} must be equal to the number of columns in \code{dmatrix}, or it can be a single character string to apply the same preference function to all criteria. Available preference functions are: "linear", "usual", "quasi-criterion", "v-shape", and "level".
  }
  \item{thr}{
    A list of threshold values for the preference functions. The structure depends on the chosen preference function for each criterion. It should be a list with the same length as the number of criteria. Each element of the list can be a named list or a numeric vector containing the required threshold parameters (e.g., \code{list(p = 50)} for "linear" and "v-shape", \code{list(s = 10)} for "quasi-criterion", \code{list(q = 5, p = 20)} for "level"). If a single list or numeric vector is provided, it will be applied to all criteria. If \code{NULL}, default thresholds based on the maximum absolute difference for each criterion will be used (with a warning).
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
  The PROMETHEE II method calculates a complete preorder of alternatives. It extends PROMETHEE I by calculating the net outranking flow (Phi_net) for each alternative, which is the difference between positive (Phi_plus) and negative (Phi_minus) outranking flows. Then, alternatives can be ranked based on their Phi_net values.
}
\value{
  A list containing the following components:
  \item{Phi_plus}{A numeric vector of positive outranking flows for each alternative.}
  \item{Phi_minus}{A numeric vector of negative outranking flows for each alternative.}
  \item{Phi_net}{A numeric vector of net outranking flows for each alternative.}
  \item{rank}{A numeric vector representing the complete ranking of alternatives based on the descending order of \code{Phi_net} values (higher is better). Ties are handled using the \code{"average"} method.}
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
  \code{\link{promethee1}}, \code{\link{promethee3}}, \code{\link{promethee4}}, \code{\link{promethee5}}, \code{\link{promethee6}}
}
\author{
 Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(
  7, 12, 8, 6, 14,   
  12, 18, 11, 8, 22, 
  8, 14, 11, 7, 10,  
  1.2, 1.9, 1.2, 0.8, 1.1, 
  0, 0, 1, 0, 1,  
  0, 1, 0, 1, 2  
), 
nrow = 5, byrow = FALSE)

colnames(dmat) <- c("C1", "C2", "C3", "C4", "C5", "C6")
rownames(dmat) <- c("A", "B", "C", "D", "E")

# Benefit-Cost vector
bc <- c(1, 1, 1, 1, -1, 1)

# Weights
uw <- c(0.13, 0.08, 0.20, 0.17, 0.12, 0.3)

# Threshold values
repeat_list <- list(list(p = 50), list(p = 30))
tvals <- rep(repeat_list, length.out = ncol(dmat))

# Preference function
pf <- rep("linear", ncol(dmat))

# Run PROMETHEE II function
respromethee2 <- promethee2(dmatrix=dmat, bcvec=bc, weights=uw,
    thr = tvals, prefuncs = pf)

# Phi+ values
print(respromethee2$Phi_plus)
 
# Phi- values
print(respromethee2$Phi_minus)

# Phi_net values
print(respromethee2$Phi_net)

# Complete ranking (PROMETHEE II)
print(respromethee2$rank)

# Run with a single threshold value applied to all criteria (linear)
tvals_single <- list(p = 40)
respromethee2_single <- promethee2(dmatrix=dmat, bcvec=bc, weights=uw,
    thr = tvals_single, prefuncs = "linear")

# Phi_net values for single threshold
print(respromethee2_single$Phi_net)

# Complete ranking for single threshold
print(respromethee2_single$rank)

# Run with different preference functions and mixed thresholds
tvals_mixed <- list(list(q = 10, p = 30), list(p = 20), list(q = 10, p = 30), 
   list(q = 10, p = 30), list(q = 10, p = 30), list(q = 10, p = 30))
pf_mixed <- c("linear", "linear", "linear", "linear", "level", "level")
respromethee2_mixed <- promethee2(dmatrix=dmat, bcvec=bc, weights=uw,
   thr = tvals_mixed, prefuncs = pf_mixed)

# Phi_net values with mixed thresholds
print(respromethee2_mixed$Phi_net)

# Complete ranking with mixed thresholds
print(respromethee2_mixed$rank)
}