\name{rankwssim}
\alias{rankwssim}
\encoding{UTF-8}
\title{Wojciech Sałabun Coefficient of Rankings Similarity (WS)}
\description{
This function calculates the Wojciech Sałabun Coefficient of Rankings Similarity (WS) between different rankings of examined MCDA methods. The WS coefficient quantifies the similarity between two rankings of MCDA methods, giving higher priority to the similarity of positions at the top of the ranking.
}
\usage{
rankwssim(rankmat)
}
\arguments{
\item{rankmat}{A numerical matrix containing the rankings. Rows represent different rankings (e.g., results from different MCDM methods), and columns represent the ranked items (e.g., alternatives). The row names of the matrix should be the names of the rankings.}
}
\details{
The \code{wranksim} function takes a rank matrix as input and computes the WS coefficient for each pair of rankings within the matrix. The WS coefficient, introduced by Sałabun and Urbaniak, is calculated using the following formula:
\eqn{WS = 1 - \sum_{i=1}^{n} \left( 2^{-R_{xi}} \cdot \frac{|R_{xi} - R_{yi}|}{\max(|1 - R_{xi}|, |N - R_{xi}|)} \right)}
where:
\describe{
\item{\eqn{WS}}{represents the value of the similarity coefficient.}
\item{\eqn{N}}{denotes the length of the ranking (number of alternatives).}
\item{\eqn{R_{xi}}}{represents the rank of the \eqn{i}-th element in ranking \eqn{x}.}
\item{\eqn{R_{yi}}}{represents the rank of the \eqn{i}-th element in ranking \eqn{y}.}
}
The WS coefficient ranges from 0 to 1, where 1 indicates perfect similarity and 0 indicates complete dissimilarity between the two rankings. This metric prioritizes the agreement in the top ranks, which is often crucial in Multiple Criteria Decision Making (MCDM) for identifying the most preferred alternatives.

Due to the dependence of the WS coefficient value on the reference dataset (\eqn{x}), the resulting matrix of WS coefficient values is asymmetrical, unlike Spearman's rank correlation coefficient. Therefore, for accurate and fair comparisons, each MCDM method's ranking set should be used as the reference (\eqn{x}) when determining the WS coefficient value.

Based on Sałabun and Urbaniak (2020), Rachman et al (2024) states that two rankings are considered to have high ranking similarity if \eqn{WS > 0.689}, while \eqn{WS < 0.352} indicates low similarity. MCDA methods are assumed to have a medium level of likeness if the WS value falls within the range from 0.352 to 0.689.
Two MCDA methods are defined having totally similar if they reach absolute consensus by achieving Spearman Rank Correlation = +1 with WS = 1.
}
\value{
The function returns a data frame containing the WS coefficients between all pairs of input rankings. The row and column names correspond to the names of the rankings from the input matrix. Each cell \((i, j)\) in the data frame represents the WS coefficient calculated with the \(i\)-th ranking as the reference (\eqn{x}) and the \(j\)-th ranking as the compared ranking (\eqn{y}). The coefficients are rounded to two decimal places and formatted as character strings.
}
\references{
Rachman, A. P., Ichwania, C., Mangkuto, R. A., Pradipta, J., Koerniawan, M. D., & Sarwono, J. (2024). Comparison of multi-criteria decision-making methods for selection of optimum passive design strategy. \emph{Energy and Buildings, 314}, 114285. <doi:10.1016/j.enbuild.2024.114285>

Sałabun, W., Urbaniak, K. (2020). A New Coefficient of Rankings Similarity in Decision-Making Problems. In: Krzhizhanovskaya, V., et al. Computational Science – ICCS 2020. ICCS 2020. Lecture Notes in Computer Science, vol 12138. Springer, Cham. <doi:10.1007/978-3-030-50417-5_47>
}
\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},\code{\link{rankentperm}},
  \code{\link{rankmia}}, \code{\link{rankwilcox}},  
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
MAUT              3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
MEGAN             5  3.0 11  7  1.0  6  2  9.0  4.0  10  8.0 12.0
MOORA             7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
PROMET1           4 10.0  6  8  1.0  5  3  9.0  2.0   7 11.0 12.0
PROMET2           3  9.0  6  8  1.0  5  4 10.0  2.0   7 11.0 12.0
PROMET6           4 10.0  6  8  1.0  5  3  9.0  2.0   7 11.0 12.0
RAM               7  1.0 11  9  2.0  4  3  8.0  6.0  10  5.0 12.0
ROV               3 12.0  4  8  1.0  7  6  9.0  2.0   5 10.0 11.0
SMART             5  3.0 11  7  2.0  6  1  9.0  4.0  10  8.0 12.0
TOPSIS            2 12.0  4  8  1.0  7  6  9.0  3.0   5 11.0 10.0
VIKOR             1 12.0  4  8  2.0  5  7 10.0  3.0   6  9.0 11.0
WASPAS            2 12.0  1  8  4.0  7  6 10.0  3.0   5 11.0  9.0
WPM               2 12.0  1  7  4.0  8  6 10.0  3.0   5 11.0  9.0
WSM               4  7.0  9  5  1.0  6  2  8.0  3.0  11 10.0 12.0
", header = TRUE, row.names = 1))

print(rankmat)

# Calculate the WS similarity matrix using rank matrix
reswsm <- rankwssim(rankmat)
print(reswsm$displaymat)
}