\name{rankspearman}
\alias{rankspearman}
\encoding{UTF-8}
\title{Spearman's Rank Correlation Matrix}
\description{
Calculates Spearman's rank correlation matrix for a ranking matrix.
}
\usage{
rankspearman(rankmat)
}
\arguments{
\item{rankmat}{
A numerical matrix containing ranking values. The function automatically transposes this matrix before calculating correlations, meaning each column should represent a criterion/variable.
}
}
\details{
In the resulting matrix, statistical significance stars appear above the diagonal, while correlation coefficients appear below the diagonal. 
Elements on the diagonal show the correlation coefficient (usually 1.00).
The function calculates Spearman's rank correlation and associated p-values for each pair of variables using \code{cor.test(method = "spearman")}.
Significance stars are assigned based on p-values:
\describe{
\item{***}{$p < 0.001$}
\item{**}{$p < 0.01$}
\item{*}{$p < 0.05$}
\item{Blank}{$p >= 0.05$ (not significant)}
}
The stars appear in cells above the diagonal show the signifiance of correlations. Cells below the diagonal contain rounded correlation coefficients. Diagonal cells show a variable's correlation with itself (1.00).
}
\value{
Returns a list:
\item{cormat}{
A formatted character matrix displaying correlation coefficients (below the diagonal) and statistical significance stars (above the diagonal).
Diagonal elements are displayed as 1.00.
}
\item{p_values}{
A numerical matrix containing calculated p-values for each variable pair.
Diagonal elements will be \code{NA} (Not Available).
}
}
\note{
If ties exist in the ranking matrix, the \code{cor.test} function may be unable to compute exact p-values and will instead use approximate values.
This usually does not affect the validity of the analysis but is important to note.
}
\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},\code{\link{rankentperm}},
  \code{\link{rankmia}}, \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankrangesim}}
}

\author{
Cagatay Cebeci
}
\examples{
# Sample ranking matrix
rankmat <- as.matrix(read.table(text = "
         G1   G2 G3 G4   G5 G6 G7   G8   G9 G10  G11  G12
ARAS       3 12.0  1  8  4.0  7  6 10.0  2.0   5 11.0  9.0
EDAS       3 12.0  2  8  4.0  9  6 10.0  1.0   5 11.0  7.0
ELECTRE4   3 11.5  4  8  1.5  7  5  9.5  1.5   6 11.5  9.5
FUCA       4  3.0 10  7  2.0  6  1  9.0  5.0  11  8.0 12.0
GRA        6  3.0 10  5  2.0  7  1  9.0  4.0  11  8.0 12.0
MABAC      3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
MAIRCA    10  1.0  9  5 12.0  6  7  4.0 11.0   8  3.0  2.0
MARCOS     3 12.0  1  8  4.0  7  6 10.0  2.0   5 11.0  9.0
MAUT       3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
MEGAN      5  3.0 11  7  1.0  6  2  9.0  4.0  10  8.0 12.0
MOORA      7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
PROMET1    4 10.0  6  8  1.0  5  3  9.0  2.0   7 11.0 12.0
PROMET2    3  9.0  6  8  1.0  5  4 10.0  2.0   7 11.0 12.0
PROMET6    4 10.0  6  8  1.0  5  3  9.0  2.0   7 11.0 12.0
RAM        7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
ROV        3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
SMART      5  3.0 11  7  2.0  6  1  9.0  4.0  10  8.0 12.0
TOPSIS     2 12.0  4  8  1.0  7  6  9.0  3.0   5 11.0 10.0
VIKOR      1 12.0  4  8  2.0  5  7 10.0  3.0   6  9.0 11.0
WASPAS     2 12.0  1  8  4.0  7  6 10.0  3.0   5 11.0  9.0
WPM        2 12.0  1  7  4.0  8  6 10.0  3.0   5 11.0  9.0
WSM        4  7.0  9  5  1.0  6  2  8.0  3.0  11 10.0 12.0",
header = TRUE, row.names = 1))

result <- rankspearman(rankmat)

print(result$displaymat)
print(result$p_values)
}
