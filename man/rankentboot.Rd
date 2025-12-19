\name{rankentboot}
\alias{rankentboot}
\title{Bootstrap-Based Rank Entropy Comparison}
\description{
  Computes pairwise entropy differences for ranking matrices using Shannon entropy or Jensen-Shannon divergence.
  The function applies a bootstrap-based statistical test and adjusts p-values using various correction methods.
}
\usage{
rankentboot(rankmat, nboot = 1000, entropyopt = "shannon", padjmethod = "none")
}
\arguments{
  \item{rankmat}{
    A ranking matrix where rows represent methods and columns represent ranking positions.
  }
  \item{nboot}{
    The number of bootstrap resampling iterations. Defaults to 1000.
  }
  \item{entropyopt}{
    The entropy calculation method. Options:
    \describe{
      \item{\code{"shannon"}}{Uses Shannon entropy.}
      \item{\code{"jsd"}}{Uses Jensen-Shannon divergence.}
    }
  }
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
  This function computes entropy-based differences between ranking methods and uses bootstrap sampling to generate p-values.
  P-values can be adjusted using multiple correction methods. The output matrix contains:
  
  \itemize{
    \item Lower triangle: Similarity scores (inverse entropy differences).
    \item Upper triangle: Bootstrap p-values (adjusted if selected).
    \item Diagonal: Empty values.
  }
}
\value{
  A list containing:
  \item{entropy}{Selected entropy method ("shannon" or "jsd").}
  \item{p_adjustment}{Selected p-value adjustment method.}
  \item{nboot}{Number of bootstrap iterations performed.}
  \item{displaymat}{Final comparison matrix with similarity scores and adjusted p-values.}
}
\seealso{
  \code{\link{p.adjust}}, \code{\link{rankcompare}}, 
  \code{\link{ranksrd}}, \code{\link{rankentperm}}, 
  \code{\link{rankmia}}, \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}, \code{\link{rankrangesim}}
}
\author{
  Cagatay Cebeci
}
\examples{

# Rank matrix
rankmat <- as.matrix(read.table(text = "
                     G1   G2 G3 G4   G5 G6 G7   G8   G9 G10  G11  G12
ARAS              3 12.0  1  8  4.0  7  6 10.0  2.0   5 11.0  9.0
EDAS              3 12.0  2  8  4.0  9  6 10.0  1.0   5 11.0  7.0
ELECTRE4          3 11.5  4  8  1.5  7  5  9.5  1.5   6 11.5  9.5
FUCA              4  3.0 10  7  2.0  6  1  9.0  5.0  11  8.0 12.0
GRA               6  3.0 10  5  2.0  7  1  9.0  4.0  11  8.0 12.0
MABAC             3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
MAIRCA           10  1.0  9  5 12.0  6  7  4.0 11.0   8  3.0  2.0
MARCOS            3 12.0  1  8  4.0  7  6 10.0  2.0   5 11.0  9.0
MEGAN             5  3.0 11  7  1.0  6  2  9.0  4.0  10  8.0 12.0
MOORA             7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
PROMET2           3  9.0  6  8  1.0  5  4 10.0  2.0   7 11.0 12.0
RAM               7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
ROV               3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
SMART             5  3.0 11  7  2.0  6  1  9.0  4.0  10  8.0 12.0
TOPSIS            2 12.0  4  8  1.0  7  6  9.0  3.0   5 11.0 10.0
VIKOR             1 12.0  4  8  2.0  5  7 10.0  3.0   6  9.0 11.0
WASPAS            2 12.0  1  8  4.0  7  6 10.0  3.0   5 11.0  9.0
", header = TRUE, row.names = 1))

print(rankmat)

res_shannon_adj_none <- rankentboot(rankmat, nboot = 100, 
   entropyopt = "shannon", padjmethod="none")
print(res_shannon_adj_none)

res_shannon_adj_bonf <- rankentboot(rankmat, nboot = 100, 
   entropyopt = "shannon", padjmethod="bonferroni")
print(res_shannon_adj_bonf)

res_jsd_adj_fdr <- rankentboot(rankmat, nboot = 100, 
   entropyopt = "jsd", padjmethod="fdr")
print(res_jsd_adj_fdr)

res_jsd_adj_holm <- rankentboot(rankmat, nboot = 100, 
   entropyopt = "jsd", padjmethod="holm")
print(res_jsd_adj_holm)
}
