\name{rankgamma}
\alias{rankgamma}
\encoding{UTF-8}
\title{Compute Goodman & Kruskal's Gamma Statistic}
\description{
  Computes the Goodman & Kruskal's Gamma statistic between all pairs of rankings in a given matrix.
}
\usage{
rankgamma(rankmat)
}
\arguments{
  \item{rankmat}{A numeric matrix where rows represent different ranking entities and columns represent ranking criteria.}
}
\value{
  A list with the following components:
  \item{rankmat}{The input ranking matrix.}
  \item{gammamat}{A matrix of Goodman & Kruskal's Gamma statistics computed between each pair of rows in the ranking matrix.}
}
\examples{
# Rank matrix based on MCDA methods
rankmat <- matrix(c(
    3, 12.0, 1, 8, 4.0, 7, 6, 10.0, 2.0, 5, 11.0, 9.0,
    3, 12.0, 2, 8, 4.0, 9, 6, 10.0, 1.0, 5, 11.0, 7.0,
    3, 11.5, 4, 8, 1.5, 7, 5, 9.5, 1.5, 6, 11.5, 9.5,
    4, 3.0, 10, 7, 2.0, 6, 1, 9.0, 5.0, 11, 8.0, 12.0,
    6, 3.0, 10, 5, 2.0, 7, 1, 9.0, 4.0, 11, 8.0, 12.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    10, 1.0, 9, 5, 12.0, 6, 7, 4.0, 11.0, 8, 3.0, 2.0,
    3, 12.0, 1, 8, 4.0, 7, 6, 10.0, 2.0, 5, 11.0, 9.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    5, 3.0, 11, 7, 1.0, 6, 2, 9.0, 4.0, 10, 8.0, 12.0,
    7, 1.0, 11, 9, 2.0, 4, 3, 8.0, 6.0, 10, 5.0, 12.0,
    4, 10.0, 6, 8, 1.0, 5, 3, 9.0, 2.0, 7, 11.0, 12.0,
    3, 9.0, 6, 8, 1.0, 5, 4, 10.0, 2.0, 7, 11.0, 12.0,
    4, 10.0, 6, 8, 1.0, 5, 3, 9.0, 2.0, 7, 11.0, 12.0,
    7, 1.0, 11, 9, 2.0, 4, 3, 8.0, 6.0, 10, 5.0, 12.0,
    3, 12.0, 4, 8, 1.0, 7, 6, 9.0, 2.0, 5, 10.0, 11.0,
    5, 3.0, 11, 7, 2.0, 6, 1, 9.0, 4.0, 10, 8.0, 12.0,
    2, 12.0, 4, 8, 1.0, 7, 6, 9.0, 3.0, 5, 11.0, 10.0,
    1, 12.0, 4, 8, 2.0, 5, 7, 10.0, 3.0, 6, 9.0, 11.0,
    2, 12.0, 1, 8, 4.0, 7, 6, 10.0, 3.0, 5, 11.0, 9.0,
    2, 12.0, 1, 7, 4.0, 8, 6, 10.0, 3.0, 5, 11.0, 9.0,
    4, 7.0, 9, 5, 1.0, 6, 2, 8.0, 3.0, 11, 10.0, 12.0
), nrow=22, byrow=TRUE) # Changed nrow to 22 as there are 22 rows of data

rownames(rankmat) <- c("ARAS", "EDAS", "ELECTRE4", "FUCA", "GRA", "MABAC", "MAIRCA", 
    "MARCOS", "MAUT", "MEGAN", "MOORA", "PROMET1", "PROMET2", "PROMET6", "RAM", "ROV",
    "SMART", "TOPSIS", "VIKOR", "WASPAS", "WPM", "WSM")
colnames(rankmat) <- c("G1", "G2", "G3", "G4", "G5", 
    "G6", "G7", "G8", "G9", "G10", "G11", "G12")

print(rankmat)

result <- rankgamma(rankmat)
print(round(result$gammamat, 2))
gammamat <- round(result$gammamat, 2)

corplot(gammamat, labsize = 10, fsize = 3, fcolor = "black", 
   colpal = c("red", "white", "dodgerblue"), title = "Goodman- Kruskall's Gamma Statistics", 
   xlab = "Methods", ylab = "Methods", shape = "square", displabels = TRUE)
}
\author{
  Cagatay Cebeci
}




