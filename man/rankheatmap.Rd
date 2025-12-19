\name{rankheatmap}
\alias{rankheatmap}
\encoding{UTF-8}
\title{Rank Heatmap Plot}
\usage{
  rankheatmap(rankmat, colpal = NULL, yl = NULL, xl = NULL, htitle = NULL, 
    tcol = "white", dendro = "row", cellnotes = FALSE)
}
\arguments{
  \item{rankmat}{A data frame or numeric matrix containing performance ranks of the alternatives (models, methods, etc.). 
    The alternatives are located in the columns of the data frame.}
  \item{colpal}{A color palette to be used in drawing the heatmap, or an integer between 1 and 6 for built-in palettes.}
  \item{xl}{A character string describing the x-axis label.}
  \item{yl}{A character string describing the y-axis label.}
  \item{htitle}{A character string describing the main title of the heatmap.}
  \item{tcol}{A character string describing the color of text in the cells of the heatmap.}
  \item{dendro}{A character string that stands for row and column dendrograms. Use 'row', 'column', or 'both'.}
  \item{cellnotes}{A logical flag indicating whether the cell notes are displayed or not. Default value is FALSE. Avoid using it for large matrices.}
}
\description{
  Draws a heatmap of the overall ranks of \code{n} number of alternatives across \code{p} metrics.
}
\details{
  A heatmap is a type of graphic where the values of elements in a matrix are represented by colors. It visualizes a data matrix with two dimensions, applying gradual colors to indicate different values of the items. Typically, darker colors signal higher values while lighter colors denote lower values. Heatmaps are popular in exploratory data analysis and data science applications. It's essential to decide on proper color encoding and meaningful reordering of rows and columns to yield informative heatmaps (Gehlenborg & Wong, 2012).

  Heatmaps are often used to identify patterns, correlations, and outliers in datasets at a glance. They find applications in fields like finance, geography, data analytics, and genomics. Additionally, they can be extremely useful for model evaluation, letting users quickly gauge the performance of alternative models based on multiple criteria.
}
\value{
 No value returned.
}
\references{
Gehlenborg, N., & Wong, B. (2012). Heat maps. \emph{Nature Methods, 9}(3), 213.
}
\seealso{
 \code{\link{boxplotmcda}}, \code{\link{rankcompare}}
}

\author{
Cagatay Cebeci
}

\examples{
# Example rank matrix
alternatives <- paste0("Alt", LETTERS[1:8])
criteria <- paste0("C", 1:15)

set.seed(24)
rankmat <- data.frame(Method = alternatives)
for (criterion in criteria) {
  rankmat[[criterion]] <- sample(1:8, 8, replace = FALSE)
}
rownames(rankmat) <- rankmat[,1]
rankmat <- as.matrix(rankmat[,-1])
print(rankmat)
rankheatmap(rankmat)

# Heatmap with color palette number 1 and row dendrogram only
rankheatmap(rankmat, colpal = 1, dendro = "row", cellnotes = TRUE)

# Heatmap with user-defined color palette, column dendrogram, and cell notes 
colpalcustom <- colorRampPalette(c("orange", "green", "blue"))
rankheatmap(rankmat, colpal = colpalcustom, dendro = "row", cellnotes = TRUE)
}
\concept{heatmaps}
\concept{model evaluation}
\concept{performance metrics}
\concept{multi-criteria decision making}
\keyword{heatmap}
\keyword{criterion}
\keyword{metrics}
\keyword{performance}
\keyword{ranking}
\keyword{plotting}
