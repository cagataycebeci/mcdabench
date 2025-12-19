\name{ranksrd}
\alias{ranksrd}
\encoding{UTF-8}
\title{SRD: Sum of Ranking Differences}
\description{
Calculates the Sum of Ranking Differences (SRD) for a set of alternatives across multiple criteria. The SRD method assesses the performance of alternatives by comparing their actual rankings to a defined
reference ranking. A lower SRD value indicates better agreement with the reference and thus, better performance.
}
\usage{
ranksrd(rankmat, refopt = "mean", defrank = NULL, tiesmethod="average")
}
\arguments{
  \item{rankmat}{
    A \code{matrix} or \code{data.frame} where rows represent alternatives and columns represent criteria. Each cell contains the rank of an alternative for a specific criterion (e.g., 1 for best, 8 for worst).
  }
  \item{refopt}{
    A character string specifying the type of reference ranking to be used. Options are:
    \itemize{
      \item \code{"mean"}: The reference rank for each alternative is its mean rank across all criteria. This is the default.
      \item \code{"best"}: The reference rank for each alternative is 1 (ideal best performance in every criterion).
      \item \code{"worst"}: The reference rank for each alternative is the maximum possible rank found within the \code{rankmat} (ideal worst performance in every criterion).
      \item \code{"custom"}: The reference rank is provided by the user via the \code{defrank} argument. This allows for a specific, externally defined ideal ranking.
    }
  }
  \item{defrank}{
    A numeric vector providing custom reference ranks for each alternative.
    This argument is required if \code{refopt} is \code{"custom"}.
    If provided, its names must match the row names of \code{rankmat}.
    Defaults to \code{NULL} when not used.
  }
  \item{tiesmethod}{The ties method to be applied in rank calculation. The default value is 'average'. Other options are 'min', 'max', and 'none'.} 
}
\details{
The SRD method, as proposed by Károly Héberger (Héberger, 2010), offers a fair and unambiguous way to compare methods, models, or alternatives. It quantifies the "goodness" of an alternative by summing the absolute 
differences between its actual rank in each criterion and its corresponding reference rank. The procedure is repeated for each alternative, and the resulting SRD values are used to rank them.

The formula for SRD for a single alternative is:
\deqn{SRD_i = \sum_{j=1}^{m} |R_{ij} - R_{ref,i}|}
where \eqn{R_{ij}} is the rank of alternative \eqn{i} in criterion \eqn{j}, \eqn{R_{ref,i}} is the reference rank for alternative \eqn{i}, and \eqn{m} is the total number of criteria.
A smaller SRD value indicates that the alternative's ranking is closer to the reference ranking, implying better performance.
}
\value{
A list containing the following components:
  \item{srdvals}{A named numeric vector of calculated SRD values for each alternative.}
  \item{srdrank}{A named numeric vector representing the ranking of alternatives based on their SRD values, from lowest (best) to highest (worst).}
  \item{refrank}{A named numeric vector of the reference ranks used for each alternative.}
}
\references{
Héberger, K. (2010). Sum of Ranking Differences (SRD) for method comparison and its relationship to the principle of parsimony. \emph{Chemometrics and Intelligent Laboratory Systems}, 101(1), 1-8.
}
\author{
Cagatay Cebeci
}
\note{
Ensure that the input \code{rankmat} contains actual rank values (e.g., 1, 2, 3...) and not raw performance scores. 
The ranks should be consistent (e.g., 1 is always best, higher numbers are worse).
}
\seealso{
\code{\link{rank}}
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

# Using mean rank as reference (default)
ressrd_mean <- ranksrd(rankmat, refopt = "mean")
print(ressrd_mean)

# Using the best possible rank (rank 1) as reference
ressrd_best <- ranksrd(rankmat, refopt = "best")
print(ressrd_best$srdvals)
print(ressrd_best$srdrank)
print(ressrd_best$rank)

# Using the worst possible rank (max rank in data) as reference
ressrd_worst <- ranksrd(rankmat, refopt = "worst")
print(ressrd_worst$srdvals)
print(ressrd_worst$srdrank)

# Using user defined reference
usrrank <- c(4, 6, 1, 8, 2, 3, 7, 5)
names(usrrank) <- rownames(rankmat) 

ressrd_user <- ranksrd(rankmat, refopt = "custom", defrank = usrrank)
print(ressrd_user$srdvals)
print(ressrd_user$srdrank)
}