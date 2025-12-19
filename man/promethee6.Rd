\name{promethee6}
\alias{promethee6}
\encoding{UTF-8}
\title{PROMETHEE VI: Preference Ranking Organization METHod for Enrichment Evaluations (Method VI)}
\description{
PROMETHEE VI method for multi-criteria decision analysis, ranks alternatives based on a balance between their central tendency and variability.
}
\usage{
promethee6 (dmatrix, bcvec, weights, normethod = NULL, prefuncs = NULL, thr = NULL,
   varmethod = "abs_sum", p = NULL, q = NULL, tiesmethod="average")
}
\arguments{
\item{dmatrix}{A numeric matrix representing the decision matrix. Rows represent alternatives, and columns represent criteria.}
\item{bcvec}{A numeric vector indicating the benefit/cost orientation of each criterion. Use 1 for benefit criteria (to be maximized) and -1 for cost criteria (to be minimized).}
\item{weights}{A vector of weights representing the relative importance of each criterion. Defaults to \code{"equal"}. See \code{\link{calcweights}} for details.}
\item{normethod}{An optional character string specifying the normalization method for the decision matrix because there is usually no need to normalize the decision matrix with PROMETHEE V. If not \code{NULL}, the \code{calcnormal} function will be used. Common options might include "linear", "vector", etc. Defaults to \code{NULL}, no normalization is performed. See \code{\link{calcnormal}} for details.}
\item{prefuncs}{Optional. A character vector specifying the preference function for each criterion. Can be a single string (applied to all criteria) or a vector with length equal to the number of criteria. Available preference functions are: "linear", "usual", "quasi-criterion", "v-shape", and "level". If \code{NULL}, "linear" is used for all criteria.}
\item{thr}{Optional. A list of threshold parameters for the preference functions. If a single list, it's applied to all criteria. If a list of lists, each element corresponds to a criterion. The structure of each element depends on the preference function:
   \itemize{
     \item{"linear", "v-shape":} A list with a single numeric value named \code{p} (preference threshold).
     \item{"quasi-criterion":} A list with a single numeric value named \code{s} (indifference threshold).
     \item{"level":} A list with two numeric values named \code{p} (preference threshold) and \code{q} (indifference threshold).
   }
   If \code{NULL}, default thresholds are calculated based on the maximum absolute difference for each criterion.}
\item{p}{A numeric value representing the preference parameter used in the calculation of Lambda.  Must be non-negative.}
\item{q}{A numeric value representing the indifference parameter used in the calculation of Lambda. Must be non-negative.}
\item{varmethod}{Variability method. The \code{"abs_sum"} for absolute sum. The \code{"range"} is also available as an option.}
\item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The \code{promethee6} function implements the PROMETHEE VI method, an extension of the PROMETHEE II method.  PROMETHEE VI considers both the central tendency and the variability of the net outranking flows to rank alternatives.

The \code{prefuncs} argument allows for different preference functions to be used for each criterion.  The preference function transforms the difference between the evaluations of two alternatives on a given criterion into a preference degree.

The \code{thr} argument specifies the parameters for the preference functions.  The required parameters vary depending on the chosen preference function.

The arguments \code{p} and \code{q} are specific to PROMETHEE VI.  They represent the preference and indifference parameters, respectively, used to calculate the Lambda values, which adjust the final ranking based on the variability of the net outranking flows.
}
\value{
A list containing the following components:
\item{Phi_plus}{A numeric vector of positive outranking flows (\eqn{\phi^+}) for each alternative.}
\item{Phi_minus}{A numeric vector of negative outranking flows (\eqn{\phi^-}) for each alternative.}
\item{Phi_net}{A numeric vector of net outranking flows (\eqn{\phi_{net}}) for each alternative.}
\item{Central_tendency}{A numeric vector representing the central tendency of the net outranking flows (equal to Phi_net).}
\item{Variability}{A numeric vector representing the variability of the outranking flows, calculated as the absolute sum of positive and negative flows.}
\item{Lambda}{A numeric vector of Lambda values for each alternative, based on its variability and the \code{p} and \code{q} parameters.}
\item{Xi}{A numeric vector of final values for each alternative, calculated by adjusting the central tendency with Lambda and variability.}
\item{rank}{A numeric vector of the final ranks of the alternatives, with lower ranks indicating better performance, based on the Xi values.}
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
  \code{\link{promethee1}}, \code{\link{promethee2}}, \code{\link{promethee3}}, \code{\link{promethee4}}, \code{\link{promethee5}}
}
\author{
  Cagatay Cebeci
}
\examples{
# Decision matrix
dmat <- matrix(c(5, 150, 3, 200, 4, 180), nrow = 3, byrow = TRUE)
rownames(dmat) <- paste0("A", 1:3)
colnames(dmat) <- paste0("C", 1:2)

# Benefit-Cost vector
bc <- c(-1, 1)

# Predefined weights
uw <- c(0.5, 0.5)

# Threshold values
tvals <- list(list(p = 50), list(p = 30))

# Preference function
pf <- c("linear", "linear")

# PROMETHEE VI p and q parameters
pval <- 120
qval <- 20

# Run PROMETHEE VI function
respromethee6_1 <- promethee6(dmatrix=dmat, bcvec=bc, weights=uw,
     prefuncs = pf, thr = tvals, p = pval, q = qval)
print(respromethee6_1)

# Using different preference functions
pf2 <- c("v-shape", "level")
thr2 <- list(list(p = 60), list(q = 15, p = 40))
respromethee6_2 <- promethee6(dmatrix=dmat, bcvec=bc, weights=uw, 
      prefuncs = pf2, thr = thr2, p = pval, q = qval)
print(respromethee6_2)

# Different weights and normalization
uw <- c(0.7, 0.3)
respromethee6_3 <- promethee6(dmatrix=dmat, bcvec=bc, weights=uw, 
  normethod = "ratio", p = pval, q = qval)
print(respromethee6_3)
}
