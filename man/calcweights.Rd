\name{calcweights}
\alias{calcweights}
\title{Calculation of Criteria Weights}
\description{
  This function calculates criteria weights using different methods such as Equal weighting, and weighting based on Angle, CRITIC, CILOS, Entropy, Geometric mean, GINI coefficient, MEREC, MPSI, 
  ROC, RS, Standard deviation, and Variance.
}
\usage{
calcweights(dmatrix, bcvec = NULL, type = "equal", normethod = NULL)
}
\arguments{
  \item{dmatrix}{A normalized decision matrix.}
  \item{bcvec}{Benefit-Cost vector, used with certain methods such as MEREC.}
  \item{type}{The type of weighting method to be applied. The options include: 
   \code{"angle"}, \code{"critic"}, \code{"cilos"}, \code{"entropy"}, \code{"equal"}, \code{"geom"}, \code{"gini"}, 
   \code{"merec"}, \code{"roc"}, \code{"rs"}, \code{"sdev"}, and \code{"var"}. Default weighting method is "equal".}
  \item{normethod}{Optinal argument for normalization method to apply. The options include:
    \code{"max"}, \code{"maxmin"}, \code{"vector"}, \code{"sum"},
    \code{"enhanced"}, \code{"linear"}, \code{"nonlinear"},
    \code{"logarithmic"}, \code{"zavadskas"}. Default is \code{NULL}. }
}
\details{
 In multi-criteria decision-making methods, normalization is often necessary when calculating weights, but it depends on the specific approach used. Normalization ensures that criteria with different scales become comparable, which is crucial for accurate decision analysis. For instance, the entropy method applies normalization to make data scale-independent. Similarly, the CRITIC method typically requires normalization since it considers the variability and correlation between criteria. Additionally, various normalization techniques, such as Min-Max normalization, Z-score normalization, and vector normalization, can be employed depending on the method's requirements. Min-Max normalization scales data within a fixed range, making comparisons straightforward. Z-score normalization adjusts data based on mean and standard deviation, which helps in standardizing different distributions. Vector normalization converts values into unit vectors, ensuring proportional representation across criteria. Selecting the most appropriate normalization technique is key to ensuring accurate and meaningful results in MCDA. See \code{link{normalize}} for available methods of normalization.

 The current version of function includes the following weighting methods:
  \describe{
    \item{\code{angle}}{Angular weighting method}
    \item{\code{critic}}{CRITIC: CRiteria Importance Through Intercriteria Correlation method}
    \item{\code{cilos}}{CLIOS: Criterion Impact LOSs method}
    \item{\code{entropy}}{Entropy: Shannon entropy method}
    \item{\code{entropy2}}{Entropy method II}
    \item{\code{equal}}{Equal weighting}
    \item{\code{geom}}{Geometric mean based weighting}
    \item{\code{gini}}{GINI coefficient based weighting}
    \item{\code{merec}}{MEREC: MEthod based on the Removal Effects of Criteria)}
    \item{\code{mpsi}}{MPSI: Modified Preference Selection Index}   
    \item{\code{roc}}{ROC: Rank Order Centroid weighting}
    \item{\code{rs}}{RS: Rank Sum weighting}
    \item{\code{sdev}}{Standard deviation-based weighting}
    \item{\code{var}}{Variance-based weighting}
 }
}
\value{
  A vector of weights corresponding to each column of the decision matrix.
}
\references{
Diakoulaki, D., Mavrotas, G., & Koumoutsos, K. (1995). Determining objective weights in multiple criteria problems using the CRITIC method.
\emph{Computers & Operations Research, 22}(7), 763-770. <doi:10.1016/0305-0548(94)00059-H>

Ghorabaee, M. K., Zavadskas, E. K., Olfat, L., & Turskis, Z. (2017). Multi-criteria inventory classification using a novel method of evaluation 
 based on ratios and comparisons (MEREC). \emph{Informatica, 28}(3), 599-616.

Gligoric, M., Gligoric, Z., Lutovac, S., Negovanovic, M., & Langovic, Z. (2022). Novel hybrid MPSI-MARA decision-making model for support system selection 
 in an underground mine. \emph{Systems, 10}(6), 248. <doi:10.3390/systems10060248>

Keshavarz-Ghorabaee, M., Amiri, M., Zavadskas, E. K., Turskis, Z., & Antucheviciene, J. (2021). Determination of objective weights using a new method 
based on the removal effects of criteria (MEREC). \emph{Symmetry, 13}(4), 525. <doi:10.3390/sym13040525>

Shuai, D., Zongzhun, Z., Yongji, W., & Lei, L. (2012, May). A new angular method to determine the objective weights. In \emph{ 2012 24th Chinese 
Control and Decision Conference (CCDC)} (pp. 3889-3892). IEEE.

Zavadskas, E. K., & Podvezko, V. (2016). Integrated determination of objective criteria weights in MCDM. 
  \emph{International Journal of Information Technology & Decision Making, 15}(2), 267-283. <doi:10.1142/S0219622016500036>

Zeleny, M. (1982). \emph{Multiple criteria decision making}. McGraw-Hill.

Jia, J., Fischer, G. W., & Dyer, J. S. (1998). Attribute weighting methods and decision quality in the presence of response error: a simulation study. 
\emph{Journal of Behavioral Decision Making, 11}(2), 85-105. <doi:10.1002/(SICI)1099-0771(199806)11:2<85::AID-BDM282>3.0.CO;2-K>
}

\author{
  Cagatay Cebeci
}
\examples{
# Load the dataset
data(egrids)
dmat <- egrids$dmat  # Decision matrix
bc <- egrids$bcvec   # Benefit-cost vector

# Calculate weights using the equal weighting method
calcweights(dmatrix=dmat, type="equal")

# Calculate weights using the entropy method
calcweights(dmatrix=dmat, type="entropy")

# Calculate weights using the entropy method with MaxMin normalization
calcweights(dmatrix=dmat, bcvec=bc, type="entropy", normethod="maxmin")

# Calculate weights using the std. deviation method
calcweights(dmatrix=dmat, type="sdev")

# Calculate weights using the MEREC method
calcweights(dmatrix=dmat, bcvec=bc, type="merec")

# Calculate weights using the MPSI method
calcweights(dmatrix=dmat, bcvec=bc, type="mpsi")

# Calculate weights using the angular method
calcweights(dmatrix=dmat, bcvec=bc, type="angle")
}

