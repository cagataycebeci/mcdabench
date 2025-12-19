\name{rankwilcox}
\alias{rankwilcox}
\encoding{UTF-8}
\title{Pairwise Wilcoxon Rank Sum Tests for Comparing Rows of a Matrix}
\description{
  Performs pairwise Wilcoxon Rank Sum tests /WSRT) to compare the rows of a given matrix.
  P-values are adjusted for multiple comparisons using the specified method.
}
\usage{
rankwilcox(rankmat, alpha = 0.05, padjmethod = "none")
}
\arguments{
  \item{rankmat}{A numeric matrix where rows represent different groups or conditions
    and columns represent paired observations or measurements.}
  \item{alpha}{The significance level for determining statistical significance (default: 0.05).}
  \item{padjmethod}{The method for adjusting p-values for multiple comparisons.
    See \code{\link[stats]{p.adjust.methods}} for available options (default: "none").}
}
\details{
  This function iterates through all unique pairs of rows in the input matrix
  and performs a Wilcoxon signed-rank test for each pair. The Wilcoxon
  signed-rank test is used to assess the significance of the difference
  between two related samples. The p-values obtained from these pairwise
  tests are then adjusted for multiple comparisons to control the family-wise
  error rate. The results are presented in a matrix format, with adjusted
  p-values in the upper triangle and test statistics in the lower triangle.
}
\value{
  A list containing the following components:
  \item{pvalmat}{A matrix of raw p-values from the pairwise Wilcoxon signed-rank tests.}
  \item{adjpvalmat}{A matrix of p-values adjusted for multiple comparisons.}
  \item{displaymat}{A data frame where the upper triangle contains the adjusted
    p-values (formatted to two decimal places) and the lower triangle contains
    the Wilcoxon test statistics (formatted to two decimal places). The diagonal
    is empty.}
  \item{alpha}{The significance level used.}
  \item{padjmet}{The method used for p-value adjustment.}
}
\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},\code{\link{rankentperm}},
  \code{\link{rankmia}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}, \code{\link{rankrangesim}}
}
\author{
Cagatay Cebeci
}

\examples{
# Rank matrix
rankmat <- matrix(c(
    3, 2, 3, 5, 6, 3, 3, 5, 5, 2, 5, 5, 2, 1, 3,
    7, 6, 7, 8, 8, 7, 7, 7, 7, 4, 7, 8, 7, 7, 7,
    1, 1, 1, 1, 2, 1, 1, 1, 1, 5, 1, 2, 1, 2, 1,
    8, 8, 8, 3, 5, 6, 8, 4, 3, 7, 3, 3, 8, 3, 3,
    2, 5, 2, 2, 1, 2, 2, 2, 2, 8, 2, 1, 3, 4, 2,
    4, 7, 5, 7, 7, 8, 4, 8, 8, 1, 8, 7, 6, 8, 8,
    6, 4, 6, 4, 3, 5, 6, 3, 4, 6, 4, 4, 5, 6, 6,
    5, 3, 4, 6, 4, 4, 5, 6, 6, 3, 6, 6, 4, 5, 5),
    nrow = 5, byrow = TRUE)
rownames(rankmat) <- paste0("MET", 1:nrow(rankmat))
colnames(rankmat) <- paste0("ALT", 1:ncol(rankmat))

res_wilcox <- rankwilcox(rankmat)
print(res_wilcox$displaymat)
}
