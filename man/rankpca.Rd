\name{rankpca}
\alias{rankpca}
\encoding{UTF-8}
\title{
Principal Component Analysis for Ranking Data with Geladi's Approach
}
\description{
Performs Principal Component Analysis (PCA) on a matrix of ranks, where rows represent alternatives and columns represent criteria. This function adapts Geladi's approach to interpreting PCA results in the 
context of multi-criteria performance, focusing on overall performance (PC1) and uniformity of behavior (PC2) of alternatives.
}
\usage{
rankpca(rankmat, biplot = TRUE, reverse = FALSE)
}
\arguments{
  \item{rankmat}{
    A \code{matrix} or \code{data.frame} where rows represent alternatives and columns represent criteria. Each cell contains the rank of an alternative for a specific criterion (e.g., 1 for best, 8 for worst). 
    **Row names must be provided** to identify alternatives.}
  \item{biplot}{If \code{TRUE} (default), a biplot showing the PCA scores for alternatives and loadings for criteria will be generated using the \code{factoextra} package. }
  \item{reverse}{If \code{FALSE} (default), the matrix is processed as original, but if  \code{TRUE} the values in rows are reversed by sorting in opposite direction}
} 
\details{
This function applies PCA to ranking data, offering insights into the underlying patterns and relationships among alternatives based on their performance across multiple criteria. It extends Geladi's concept of using 
PCA to analyze "regression vectors" to ranking vectors.

The interpretation of the first two principal components (PC1 and PC2) follows Geladi's framework:
\itemize{
  \item \strong{PC1 (Overall Performance):} This component typically captures the largest amount of variance and differentiates alternatives based on their overall performance. Alternatives clustering at one end of 
  PC1 generally perform better or worse across most criteria. The exact interpretation (e.g., positive PC1 score = better/worse) depends on the signs of the criteria loadings.
  \item \strong{PC2 (Uniformity of Behavior):} This component represents the consistency of an alternative's performance across different criteria. Alternatives with PC2 scores close to zero exhibit uniform behavior, 
  meaning their ranking profile is consistent. Alternatives with PC2 scores significantly deviating from zero show non-uniform behavior, indicating their performance varies considerably depending on the criterion.
}

\strong{Note:} PCA on ranking data is primarily a \strong{discovery and understanding tool}, not a direct ranking method. It helps identify clusters, outliers, and key driving criteria, complementing explicit ranking methods like SRD.
}
\value{
A \code{list} containing the following components:
  \item{pca_results}{The raw \code{prcomp} object returned by the PCA analysis.}
  \item{alternative_scores}{A \code{data.frame} of principal component scores for each alternative (rows are alternatives, columns are PCs).}
  \item{criteria_loadings}{A \code{data.frame} of loadings for each criterion on the principal components (rows are criteria, columns are PCs).}
  \item{explained_variance}{A numeric vector indicating the proportion of variance explained by each principal component.}
  \item{cumulative_variance}{A numeric vector indicating the cumulative proportion of variance explained by the principal components.}
  \item{pc1_interpretation}{A character string providing an interpretation of PC1's meaning based on the data's characteristics.}
  \item{pc2_interpretation}{A character string providing an interpretation of PC2's meaning.}
}
\references{
Geladi, P. (1995). The regression model comparison plot (REMOCOP), Proc. Spectroscopy Across the Spectrum IV,
 Norwich, UK, 11–14 July 1994, in: D. Andrews, A. Davies (Editors), Frontiers inAnalytical Spectroscopy, 
 The Royal Society of Chemistry, Cambridge, UK, 1995, pp. 225–236.

Héberger, K. (2010). Sum of Ranking Differences (SRD) for method comparison and its relationship to the principle of parsimony.
\emph{Chemometrics and Intelligent Laboratory Systems}, 101(1), 1-8.
}
\seealso{
\code{\link{prcomp}}, \code{\link[factoextra]{fviz_pca_biplot}}, \code{\link{ranksrd}}
}
\examples{
# Sample rank matrix (8 alternatives, 15 criteria)
rankmat <- data.frame(
  C1 = c(7, 3, 2, 8, 1, 6, 5, 4), C2 = c(2, 1, 4, 7, 6, 3, 8, 5),
  C3 = c(5, 8, 1, 7, 4, 6, 2, 3), C4 = c(8, 2, 1, 4, 3, 5, 6, 7),
  C5 = c(1, 5, 4, 6, 2, 3, 7, 8), C6 = c(6, 8, 2, 3, 1, 4, 5, 7),
  C7 = c(6, 7, 1, 8, 2, 3, 5, 4), C8 = c(6, 3, 4, 5, 2, 1, 7, 8),
  C9 = c(1, 8, 3, 4, 6, 5, 7, 2), C10 = c(7, 8, 5, 2, 4, 6, 3, 1),
  C11 = c(6, 2, 3, 8, 5, 7, 4, 1), C12 = c(6, 8, 3, 7, 4, 1, 5, 2),
  C13 = c(2, 3, 7, 5, 6, 1, 8, 4), C14 = c(2, 1, 4, 5, 3, 8, 6, 7),
  C15 = c(1, 6, 4, 5, 8, 7, 3, 2),
  row.names = c("AltA", "AltB", "AltC", "AltD", "AltE", "AltF", "AltG", "AltH")
)

# Perform Geladi PCA and show the biplot
respca <- rankpca(rankmat, biplot = TRUE, reverse=TRUE)

# Geladi PCA Analysis Results
print(paste("Proportion of Explained Variance:", paste(round(respca$expvar, 3), collapse = ", ")))

# Alternative Scores (PC1 and PC2)
print(round(respca$ascores[, 1:2], 2))

# Criteria Loadings (PC1 and PC2)
print(round(respca$critloads[, 1:2], 2))

# PC1 Interpretation
print(respca$intprtpc1)

# PC2 Interpretation
print(respca$intprtpc2)
}