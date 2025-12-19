\name{flowplot}
\alias{flowplot}
\encoding{UTF-8}
\title{Visualize Dominance Flow of Ranked Alternatives}
\description{
Creates a directional plot illustrating the dominance flow among alternatives based on their ranks. Equal ranks are marked with "=" and dominance with ">".
}
\usage{
flowplot(rankvec, colpal, txtcol = "black", orientation = "horizontal")
}
\arguments{
  \item{rankvec}{A numeric vector of rankings for each alternative. Names are used as labels; if missing, default labels \code{A1}, \code{A2}, ... are generated.}
  \item{colpal}{A vector of colors corresponding to unique rank values. Should have at least as many colors as unique entries in \code{rankvec}.}
  \item{txtcol}{Color for vertex labels. Defaults to black.}
  \item{orientation}{Layout direction of the flow. Use \code{"horizontal"} (default) or \code{"vertical"}.}
}
\details{
The function draws a graph in which each node represents an alternative. Arrows indicate flow of preference: a ">" means dominance, while "=" indicates tie. 
The layout is determined by the specified orientation. Edge labels are automatically positioned slightly above the edge direction based on flow geometry.
}

\examples{
rankvec <- c(1, 2, 3, 4, 5.5, 5.5, 7)
colpal <- rainbow(length(unique(rankvec)))

# Horizontal flow
flowplot(rankvec, colpal, orientation = "horizontal")

# Vertical flow
flowplot(rankvec, colpal, orientation = "vertical")
}

\author{
Cagatay Cebeci
}
