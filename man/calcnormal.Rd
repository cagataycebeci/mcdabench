\name{calcnormal}
\alias{calcnormal}
\encoding{UTF-8}
\title{Normalization of Decision Matrices}
\description{
  This function provides various normalization techniques to preprocess decision matrices for multi-criteria decision-making methods.
}
\usage{
calcnormal(X, bcvec, type = "maxmin")
}
\arguments{
  \item{X}{A decision matrix with dimensions \code{n x m} (n rows and m columns).}
  \item{bcvec}{A benefit-cost vector, where \code{-1} indicates cost criteria and \code{1} indicates benefit criteria.}
  \item{type}{The normalization method to apply. The options include:
    \code{"maxmin"}, \code{"vector"}, \code{"ratio"}, \code{"sum"},
    \code{"enhanced"}, \code{"linear"}, \code{"nonlinear"},
    \code{"logarithmic"}, \code{"zavadskas"}, \code{"null"}.}
}
\details{
  The supported normalization methods are:
  \describe{
    \item{\code{"null"}}{No normalization.}
    \item{\code{"maxmin"}}{Max-Min normalization. Based on Ersoy (2021); Baydaş et al. (2024).}
    \item{\code{"vector"}}{Vector normalization. Referenced in Nijkamp & van Delft (1977); Milani et al. (2005).}
    \item{\code{"sum"}}{Sum normalization. Jee & Kang (2000); Milani et al. (2005); ÖWang & Luo (2010); Stanujkic et al. (2013); Ersoy (2021).}
    \item{\code{"zavadskas"}}{Zavadskas-Turskis normalization. Zavadskas & Turskis (2011).}
    \item{\code{"enhanced"}}{Enhanced accuracy normalization. Proposed by Zeng et al. (2013).}
    \item{\code{"linear"}}{Linear normalization.}
    \item{\code{"nonlinear"}}{Non-linear normalization.}
    \item{\code{"logarithmic"}}{Logarithmic normalization.}
    \item{\code{"ratio"}}{Linear ratio normalization.}
  }
}
\value{
  A normalized decision matrix of the same dimension as \code{X}, the original decision matrix.
}
\references{
  Baydas, M., Yilmaz, M., Jovic, Z., Stevic, Z., Ozuyar, S. E. G., & Ozcil, A. (2024). 
  A comprehensive MCDM assessment for economic data: success analysis of maximum normalization, CODAS, and fuzzy approaches.
  \emph{Financial Innovation, 10}(1), 105.

  Ersoy, N. (2021). Selecting the best normalization technique for ROV method: Towards a real life application.
  \emph{Gazi University Journal of Science, 34}(2), 592-609.

  Jee, D. H., & Kang, K. J. (2000). A method for optimal material selection aided with decision making theory.
  \emph{Materials & Design, 21}(3), 199-206.

  Milani, A. S., Shanian, A., Madoliat, R., & Nemes, J. A. (2005). The effect of normalization norms in multiple attribute 
  decision making models: a case study in gear material selection.
  \emph{Structural and Multidisciplinary Optimization, 29}, 312-318.

  Nijkamp, P., & van Delft, A. (1977). Multi-criteria analysis and regional decision-making (Vol. 8). Springer Science & Business Media.

  Stanujkic, D., Djordjevic, B., & Djordjevic, M. (2013). Comparative analysis of some prominent MCDM methods:
  A case of ranking Serbian banks. \emph{Serbian Journal of Management, 8}(2), 213-241.

  Wang, Y. M., & Luo, Y. (2010). Integration of correlations with standard deviations for determining attribute 
  weights in multiple attribute decision making. \emph{Mathematical and Computer Modelling, 51}(1-2), 1-12.

  Zavadskas, E. K., & Turskis, Z. (2011). Multiple criteria decision making (MCDM) methods in economics: an overview.
  \emph{Technological and Economic Development of Economy, 17}(2), 397-427.

  Zeng, Q. L., Li, D. D., & Yang, Y. B. (2013). VIKOR method with enhanced accuracy for multiple criteria decision 
  making in healthcare management. \emph{Journal of Medical Systems, 37}(2), 9908.
}
\seealso{
  \code{\link{calcweights}}
}
\author{
  Cagatay Cebeci
}
\examples{
data(egrids)
dmatrix <- egrids$dmat
bcvec <- egrids$bcvec

# Max-min normalization
normmaxmin <- calcnormal(X=dmatrix, bcvec=bcvec, type="maxmin")
print(normmaxmin)

# Vector normalization
normvector <- calcnormal(X=dmatrix, bcvec=bcvec, type="vector")
print(normvector)

# Zavadskas normalization
normzavad <- calcnormal(X=dmatrix, bcvec=bcvec, type="zavadskas")
print(normzavad)
}
