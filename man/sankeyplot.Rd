\name{sankeyplot}
\alias{sankeyplot}
\encoding{UTF-8}
\title{Sankey Diagram}
\description{
  Generates a Sankey diagram to visualize the relationship between ranking methods and alternatives.
}
\usage{
sankeyplot(rankmat, fs = 12, nw = 30)
}
\arguments{
  \item{rankmat}{A numeric matrix of rankings, with rows representing ranking methods and columns representing alternatives. Lower rank indicates higher preference.}
  \item{fs}{Font size for node labels. Default is 12.}
  \item{nw}{Width of each node in pixels. Default is 30.}
}
\details{
  Ideal for MCDM contexts where ranking comparison is crucial.
  The Sankey diagram shows the mapping between ranking \emph{methods} (sources) and \emph{alternatives} (targets), where link thickness reflects the strength of preference. Each link's width is computed as \code{k - Rank}, meaning that lower ranks (i.e., better positions) result in thicker flows.

  This visualization allows quick identification of which alternatives are consistently favored across multiple methods, and highlights divergence in method behavior. For example, an alternative consistently ranked highly will attract wider links from several methods, creating a visually dominant path toward that target.

}
\value{
  An interactive Sankey diagram rendered using the \code{networkD3} package.
}
\seealso{
  \code{\link[networkD3]{sankeyNetwork}}
}
\examples{
rankmat <- matrix(c(
  5,9,7,8,4,2,10,11,3,1,12,6,
  3,12,8,4,9,6,7,10,2,5,11,1,
  2,12,4,8,6,7,9,10,5,1,11,3,
  2,11,6,8,3,5,9,10,1,4,12,7,
  7,9,6,8,5,2,11,10,3,1,12,4,
  5,9,7,8,4,2,10,11,3,1,12,6,
  5,8,6,9,4,2,10,11,3,1,12,7,
  4,6,9,8,1,2,7,10,3,5,11,12,
  5,8,6,9,1,3,7,12,2,4,11,10,
  5,7,6,9,1,2,8,12,3,4,11,10
), nrow = 10, byrow = TRUE)
rownames(rankmat) <- c("Equal","Critic","StdDev","Gini","Geometric",
  "Merec","Mpsi","Entropy","Roc","Rs")
colnames(rankmat) <- paste0("G",1:12)
sankeyplot(rankmat)
}
