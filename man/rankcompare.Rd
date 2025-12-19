\name{rankcompare}
\alias{rankcompare}
\encoding{UTF-8}
\title{Comparison of Ranking Results}
\description{
Computes and returns multiple ranking comparison measures, including Spearman correlation, Bootstrap tests based-on Jensen-Shannon divergences, Permutation test based-on Shannon-entropy differences,
Wojciech Salabun similarity coefficient, and Wilcoxon rank sum test statistics.
}
\usage{
rankcompare(rankmat, nperms = 1000, nboot = 1000, 
    entropyopt = "jsd", j = 1, k = 3, 
    alpha = 0.05, padjmethod = "fdr",
    biplot = TRUE)
}
   
\arguments{
  \item{rankmat}{A matrix containing ranking data where columns represent different ranking scenarios.}
  \item{nperms}{The number of permutations for entropy-based comparison (default: 1000).}
  \item{nboot}{The number of resampling for entropy-based comparison using bootstrap (default: 1000).}
  \item{entropyopt}{Entropy method: "shannon" for Shannon entropy or "jsd" for Jensen-Shannon divergence (default: "jsd").}
  \item{j}{For modified index of agreement, an integer specifying the exponent for the absolute difference, typically set to 1.}
  \item{k}{For rank similarity for a ranking range from k1 to k2. Usually set to 3 for top 3 rank similarity.}
  \item{alpha}{Significance threshold for statistical test such as Kruskal-Wallis test (default: 0.05).}
  \item{padjmethod}{P-value adjustment method for Kruskal-Wallis test (default: "fdr").}
  \item{biplot}{Logical whether Geladi's PCA-bilot is displayed or not.}
}
\value{
  A list with five elements:
  \item{src}{A matrix of Spearman correlation comparisons}
  \item{entperm}{A matrix of entropy-based comparisons with permutation tests}
  \item{entboot}{A matrix of entropy-based comparisons with bootstrap}
  \item{wsrs}{A matrix of Wojciech Salabun similarity comparisons}
  \item{mia}{A matrix of modified index of agreement comparisons}  
  \item{wilcox}{A matrix of Wilcoxon Rank Sum test comparisons}
  \item{rangesim}{A matrix of rank range similarities}
  \item{geladipca}{Geladi's PCA analysis results and PCA biplot.}  
}
\seealso{
  \code{\link{sensana}}, \code{\link{rankentboot}}, \code{\link{rankentperm}}, 
  \code{\link{rankmia}}, \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}, \code{\link{rankrangesim}}
}
\author{
Cagatay Cebeci
}
\examples{
rankmat <- matrix(sample(1:10, 50, replace=TRUE), ncol=5)
rownames(rankmat) <- paste0("A", 1:10)
colnames(rankmat) <- paste0("M", 1:5)

result <- rankcompare(rankmat, entropyopt="jsd", nboot=100, nperms=100, 
   padjmethod="bonferroni", k=3, j=1)
print(result)
}
