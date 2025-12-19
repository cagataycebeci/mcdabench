\name{promethee5}
\alias{promethee5}
\encoding{UTF-8}
\title{PROMETHEE V: Preference Ranking Organization METHod for Enrichment Evaluations (Method V)}
\description{
PROMETHEE V method applies multi-criteria decision analysis, which ranks alternatives considering both their net outranking flow and feasibility constraints.
}
\usage{
promethee5(dmatrix, bcvec, weights, normethod = NULL,
    prefuncs = NULL, thr = NULL, g = 0, l = 100, tiesmethod="average")
}
\arguments{
\item{dmatrix}{A numeric matrix representing the decision matrix. Rows represent alternatives, and columns represent criteria.}
\item{bcvec}{A numeric vector indicating the benefit/cost orientation of each criterion. Use 1 for benefit criteria (to be maximized) and -1 for cost criteria (to be minimized).}
\item{weights}{A vector of weights representing the relative importance of each criterion. See \code{\link{calcweights}} for details.}
\item{normethod}{An optional character string specifying the normalization method for the decision matrix because there is usually no need to normalize the decision matrix with PROMETHEE V. If not \code{NULL}, the \code{calcnormal} function will be used. Common options might include "linear", "vector", etc. Defaults to \code{NULL}, no normalization is performed. See \code{\link{calcnormal}} for details.}
\item{prefuncs}{
    A character vector specifying the preference function to be used for each criterion. The length of \code{prefuncs} must be equal to the number of columns in \code{dmatrix}, or it can be a single character string to apply the same preference function to all criteria. Available preference functions are: "linear", "usual", "quasi-criterion", "v-shape", and "level".}
\item{thr}{
    A list of threshold values for the preference functions. The structure depends on the chosen preference function for each criterion. It should be a list with the same length as the number of criteria. Each element of the list can be a named list or a numeric vector containing the required threshold parameters (e.g., \code{list(p = 50)} for "linear" and "v-shape", \code{list(s = 10)} for "quasi-criterion", \code{list(q = 5, p = 20)} for "level"). If a single list or numeric vector is provided, it will be applied to all criteria. If \code{NULL}, default thresholds based on the maximum absolute difference for each criterion will be used (with a warning).}
\item{g}{A numeric value representing the lower limit for the net outranking flow (\eqn{\phi_{net}}) for an alternative to be considered feasible. Must be non-negative. Default is 0.}
\item{l}{A numeric value representing the upper limit for the net outranking flow (\eqn{\phi_{net}}) for an alternative to be considered feasible. Must be non-negative. Default 100.}
\item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
PROMETHEE V method, an extension of the PROMETHEE II method is used for ranking alternatives in multi-criteria decision problems. PROMETHEE V incorporates feasibility constraints, defined by lower (\code{g}) and upper (\code{l}) limits on the net outranking flow.

The \code{prefuncs} argument allows for different preference functions to be used for each criterion.  The preference function transforms the difference between the evaluations of two alternatives on a given criterion into a preference degree.

The \code{thr} argument specifies the threshold parameters for the preference functions.  The required parameters vary depending on the chosen preference function.
}
\value{
A list containing the following components:
\item{Phi_plus}{A numeric vector of positive outranking flows (\eqn{\phi^+}) for each alternative.}
\item{Phi_minus}{A numeric vector of negative outranking flows (\eqn{\phi^-}) for each alternative.}
\item{Phi_net}{A numeric vector of net outranking flows (\eqn{\phi_{net}}) for each alternative.}
\item{Feasibility}{A logical vector indicating whether each alternative is feasible (TRUE) or not (FALSE) based on the \code{g} and \code{l} limits.}
\item{rank}{A numeric vector of the final ranks of the alternatives, with lower ranks indicating better performance. Infeasible alternatives are ranked lower than feasible ones.}
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
  \code{\link{promethee1}}, \code{\link{promethee2}}, \code{\link{promethee3}}, \code{\link{promethee4}}, \code{\link{promethee6}}
}
\author{
 Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(5, 150, 3, 200, 4, 180, 6, 170), nrow = 4, byrow = TRUE)
rownames(dmat) <- c("A1", "A2", "A3", "A4")
colnames(dmat) <- c("C1", "C2")

# Benefit-Cost vector
bc <- c(-1, 1)

# Weights
uw <- c(0.6, 0.4)

# Threshold values
tvals <- list(list(p = 50), list(p = 30))

# Parameters for PROMETHEE V
gval <- 0
lval <- 100

# Preference function
pf <- c("linear", "linear")

# Run PROMETHEE V function
respromethee5 <- promethee5(dmatrix=dmat, bcvec=bc, weights=uw, 
     prefunc = pf, thr = tvals, g = gval, l = lval)
print(respromethee5)
}
