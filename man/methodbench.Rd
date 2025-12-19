\name{methodbench}
\alias{methodbench}
\encoding{UTF-8}
\title{Benchmark Multiple MCDM Methods}
\description{
  Runs and compares multiple Multi-Criteria Decision-Making (MCDM) methods on a given decision matrix, using consistent 
  criteria weights and types. Supports both predefined and user-supplied methods with optional parameter customization.
}
\usage{
methodbench(dmatrix,bcvec, weights, normethod = "maxmin",
  mcdm = "all", params = NULL, tiesmethod = "average")
}
\arguments{
  \item{dmatrix}{A numeric decision matrix, where rows represent alternatives and columns represent criteria.}
  \item{bcvec}{A vector indicating whether each criterion is a benefit (\code{1}) or cost (\code{0}) type.}
  \item{weights}{A numeric vector of weights corresponding to each criterion.}
  \item{normethod}{Normalization method used by applicable methods (e.g., \code{"maxmin"}, \code{"ratio"}).}
  \item{mcdm}{A character vector of method names to execute. Use \code{"all"} to include all built-in methods. 
    User-defined functions are also supported.}
  \item{params}{An optional named list of method-specific parameter lists. Each element name should match a method name, 
    and its value is a list of arguments passed to that method.}
  \item{tiesmethod}{How ties in ranking should be resolved. Passed to \code{rank()}, e.g., \code{"average"}, 
    \code{"first"}, etc.}
}
\value{
  A list containing:
  \item{dmatrix}{The original decision matrix.}
  \item{weights}{The input weight vector.}
  \item{rankmat}{A matrix of rankings returned by each method. Rows represent methods, columns represent alternatives.}
}
\details{
  This function is designed to support benchmarking of multiple MCDM methods under consistent problem setups. 
  It supports built-in methods "aras","aroman","cocoso","codas","copras","edas","electre4", "fuca","gra","mabac",
    "macont","marcos","mairca","maut","mavt", "megan","megan2","moora","ocra","oretes", "promethee1","promethee2",
    "promethee3","promethee4","promethee5","promethee6", "ram","rov","smart","topsis","vikor","waspas","wpm","wsm", 
  and also allows integration of external functions if they return a named list including a \code{rank} vector.

  For each method, base arguments (\code{dmatrix}, \code{bcvec}, \code{weights}) are passed, and if any \code{params} are 
  specified, they are merged in before calling the method via \code{do.call()}.

  If a method fails or does not return a \code{rank} element, a warning is displayed and the result is skipped.
}
\seealso{
  \code{\link{calcweights}}, \code{\link{rankaggregate}}, \code{\link{rankcompare}}, \code{\link{rankheatmap}}
}

\author{
Cagatay Cebeci
}

\examples{
# Sample decision matrix
dm <- matrix(c(
  10, 20, 30, 1.5, 102, 55,
  15, 25, 35, 1.6, 90, 60,
  12, 22, 32, 1.7, 100, 58,
  13, 24, 33, 1.8, 95, 57,
  14, 26, 37, 1.9, 98, 59,
  11, 23, 31, 1.65, 101, 56,
  16, 27, 36, 1.55, 97, 61,
  17, 28, 38, 1.7, 99, 63,
  18, 29, 39, 1.8, 94, 62,
  19, 30, 40, 1.75, 96, 64
), nrow = 10, byrow = TRUE)
colnames(dm) <- paste0("C", 1:ncol(dm))
rownames(dm) <- paste0("ALT", 1:nrow(dm))

# Benefit-Cost vector
bc <- c(1, -1, 1, -1, 1, 1)

# User-defined weights
userwei <- c(0.3, 0.1, 0.2, 0.1, 0.2, 0.1)

prmlist <- list(
      aras      = list(),
      aroman    = list(lambda = 0.5, beta = 0.5),
      cocoso    = list(lambda = 0.5),
      codas     = list(thr = 0.1),
      copras    = list(),
      electre1  = list(),
      electre2  = list(),
      electre3  = list(),
      electre4  = list(),
      fuca      = list(),
      gra       = list(idesol = NULL, grdmethod = "sum", rho = 0.5),
      mabac     = list(),
      macont6   = list(p = 0.5, q = 0.5, delta = 0.5, theta = 0.5),
      marcos    = list(),
      mairca    = list(),
      maut      = list(utilfuncs = NULL, normutil = TRUE, ss = 1),
      mavt      = list(valfuncs = NULL, normvals = TRUE, ss = 1),
      megan     = list(normethod = "ratio", thr = 0, tht = "sdev"),
      megan2    = list(normethod = "ratio", thr = 0, tht = "sdev"),
      moora     = list(),
      ocra      = list(),
      oretes    = list(domplot = FALSE),
      promethee1 = list(),
      promethee2 = list(),
      promethee3 = list(strict = FALSE),
      promethee4 = list(alpha = 0.2),
      promethee5 = list(g = 0, l = 100),
      promethee6 = list(varmethod = "abs_sum"),
      ram       = list(normethod = "sum"),
      rov       = list(normethod = "maxmin"),
      smart     = list(),
      topsis    = list(normethod = "maxmin"),
      vikor     = list(normethod = "maxmin", v = 0.8),
      waspas    = list(normethod = "linear", v = 0.5),
      wpm       = list(normethod = "vector"),
      wsm       = list()
)

# prmlist$electre4 <- list(
#      p = c(0.6, 0.7, 0.5, 0.6, 0.6, 0.8),
#      q = c(0.3, 0.35, 0.25, 0.3, 0.3, 0.4),
#      v = c(0.9, 1.0, 0.8, 0.9, 0.9, 1.1))

# Compare all methods with default settings using user-defined weights
rescomp_1 <- methodbench(dmatrix=dm, bcvec=bc, weights=userwei, params=prmlist, mcdm="all")

print(rescomp_1$rankmat)

rankheatmap(rescomp_1$rankmat, colpal=1, cellnotes=TRUE, tcol="black")

# Overall ranks and flow of dominance
resoverall <- rankaggregate(rescomp_1$rankmat, tiesmethod="average", 
   topk = 3, damp = 0.5, niters = 200, tol = 1e-4)
str(resoverall)

# Compare selected methods with 'equal' weights
methodlist <- c("aras", "cocoso", "copras", "edas", "electre4", "mairca", 
   "mabac", "marcos", "megan", "ocra", "oretes", "promet1", "promet6", 
   "topsis", "vikor", "waspas")

equwei <- calcweights(dm, bcvec=bc, type="equal")
rescomp_2 <- methodbench(dmatrix = dm, bcvec = bc, weights = equwei,
   mcdm = methodlist, params=prmlist)

print(rescomp_2$rankmat)
rankheatmap(rescomp_2$rankmat, colpal=4, cellnotes=TRUE, tcol="white")

# Overall ranks and outranking
resoverall <- rankaggregate(rescomp_2$rankmat, tiesmethod="average", 
   topk = 3, damp = 0.5, niters = 200, tol = 1e-4)
print(resoverall$preference_ranking)
print(resoverall$preference_table)

# Compare selected methods with 'gini' weights
giniwei <- calcweights(dm, bcvec=bc, type="gini")
rescomp_3 <- methodbench(dmatrix = dm, bcvec = bc, weights = giniwei,
   mcdm = methodlist, params=prmlist)
print(rescomp_3$rankmat)
rankheatmap(rescomp_3$rankmat, colpal=1, cellnotes=TRUE, tcol="white")

# Overall ranks and outranking
resoverall <- rankaggregate(rescomp_3$rankmat, tiesmethod="average", topk = 3,
   damp = 0.5, niters = 200, tol = 1e-4)
print(resoverall)

# Custom function example
# A Dummy user-defined method
johndoe <- function(dmatrix, bcvec, weights, janefactor = 1) {
  score <- rowSums(scale(dmatrix) * weights) * janefactor
  ranking <- rank(-score, ties.method="average")
  results <- list(score = score, rank=ranking)
  return(results)
}

prmlist <- list(
   aras = list(),
   vikor = list(v = 0.5),
   waspas = list(v = 0.5),
   johndoe = list(jane = 0.25)
)

# Compare user-defined 'johndoe' with popular MCDM methods
rescomp_user <- methodbench(
  dmatrix = dm,
  bcvec = bc,
  weights = userwei,
  mcdm = c("edas", "topsis", "vikor", "aras", "waspas", "oreste", "johndoe"),
  params = prmlist
)

print(rescomp_user$rankmat)

rankheatmap(rescomp_user$rankmat, colpal = 3, cellnotes = TRUE, tcol = "black")

# Aggregate ranking
resoverall_user <- rankaggregate(rescomp_user$rankmat, tiesmethod = "average", topk = 3)
print(resoverall_user$preference_ranking)
print(resoverall_user$preference_table)

}

\keyword{MCDM methods}
\keyword{benchmarking}
\keyword{decision-making}
