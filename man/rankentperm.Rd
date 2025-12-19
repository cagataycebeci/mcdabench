\name{rankentperm}
\alias{rankentperm}
\encoding{UTF-8}
\title{Entropy-Based Pairwise Rank Comparison}
\description{
Compares pairwise rankings using either Shannon Entropy or Jensen-Shannon Divergence (JSD), with permutation testing for statistical significance.
}

\usage{
rankentperm(rankmat, nperms = 1000, entropyopt = "shannon", padjmethod="none")
}

\arguments{
  \item{rankmat}{A matrix representing ranked alternatives for multiple scenarios. Rows correspond to alternatives, and columns to scenarios.}
  \item{nperms}{An integer specifying the number of permutations for statistical testing. Default is \code{1000}.}
  \item{entropyopt}{A character string specifying the entropy method: \code{"shannon"} for Shannon Entropy or \code{"jsd"} for Jensen-Shannon Divergence. Default is \code{"shannon"}.}
  \item{padjmethod}{
    The p-value adjustment method. Options include:
    \describe{
      \item{\code{"holm"}}{Holm correction (step-down method).}
      \item{\code{"hochberg"}}{Hochberg correction (step-up method).}
      \item{\code{"hommel"}}{Hommel correction (enhanced Holm method).}
      \item{\code{"bonferroni"}}{Bonferroni correction (dividing p-values by the number of tests).}
      \item{\code{"BH"}}{Benjamini-Hochberg procedure (controls False Discovery Rate).}
      \item{\code{"BY"}}{Benjamini-Yekutieli procedure (more conservative FDR control).}
      \item{\code{"fdr"}}{False Discovery Rate correction (alias for BH method).}
      \item{\code{"none"}}{No adjustment applied (raw bootstrap p-values).}
    }
  }
}

\details{
This function analyzes ranking stability by computing entropy-based differences between rank distributions across weighting scenarios or MCDA methods.
The function allows the choice of Shannon Entropy or Jensen-Shannon Divergence, offering different sensitivity levels. Shannon Entropy was introduced the concept of entropy 
in Shannon's paper, 'A Mathematical Theory of Communication' in 1949, laying the foundation for information theory. The Jensen-Shannon Divergence (JSD) is based on the work 
of Johan Jensen and Claude Shannon. It is a symmetric version of the Kullback-Leibler Divergence: The KL Divergence is an asymmetric metric \eqn{D(P||Q) \neq D(Q||P)}  
The JSD was designed to address this asymmetry issue, making it a symmetric measure. 
A permutation test evaluates the statistical significance of obtained entropy based ranking differences.
}
\value{
A list containing:
\item{nperms}{The number of permutations used in the analysis.}
\item{entropyopt}{The entropy method used (\code{"shannon"} or \code{"jsd"}).}
\item{entropy_matrix}{A matrix storing entropy-based differences between scenarios.}
\item{pval_matrix}{A matrix containing p-values from permutation tests, indicating whether ranking distributions are significantly different.}
}
\references{
Shannon, C.E., Weaver, W. (1949). The Mathematical Theory of Communication, Univ of Illinois Press. <ISBN:0-252-72548-4>

Maasoumi, E., & Racine, J. (2002). Entropy and predictability of stock market returns. \emph{Journal of Econometrics, 107}(1-2), 291-312.
}
\note{
Please note that JSD can provide a slightly more sensitive difference assessment compared to Shannon Entropy.
}

\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},
  \code{\link{rankmia}}, \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}, \code{\link{rankrangesim}}
}

\author{
Cagatay Cebeci
}

\examples{
# Example rank matrix
rankmat <- matrix(c(
  1, 1, 1, 1, 2, 2, 1, 1, 2, 2,
  2, 2, 2, 3, 3, 1, 2, 2, 3, 3,
  3, 3, 3, 2, 1, 3, 3, 3, 1, 1,
  4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5
), nrow = 5, byrow = TRUE)

rownames(rankmat) <- c("OptA", "OptB", "OptC", "OptD", "OptE")
colnames(rankmat) <- paste0("SW_", 1:10)

# Shannon Entropy method
res_shannon <- rankentperm(rankmat, nperms=100, entropyopt = "shannon", padjmethod="fdr")
print(res_shannon$displaymat)
print(res_shannon$pval_matrix)

# Jensen-Shannon Divergence method
res_jsd <- rankentperm(rankmat, nperms=100, entropyopt = "jsd", padjmethod="bonferroni")
print(res_jsd$displaymat)
print(res_jsd$pval_matrix)
}

