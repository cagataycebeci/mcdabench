\name{genrandwei}
\alias{genrandwei}
\encoding{UTF-8}
\title{Random Modification of Initial Weights}
\description{
  Generates a set of weights by random modifying the initial weights over multiple iterations using a defined step size.
}
\usage{
genrandwei(initwei, ss = NULL, niters = 5)
}
\arguments{
  \item{initwei}{A numeric vector representing the initial weights assigned to each criterion.}
  \item{ss}{The step size used for modifying weights. If \code{NULL}, it defaults to 5\% of the smallest initial weight.}
  \item{niters}{Number of iterations to generate random weight modifications. Defaults to 5.}
}
\details{
  This function iteratively adjusts the initial weight values by adding random perturbations within the defined step size. 
  The modified weights are normalized to ensure their sum remains equal to 1. Additionally, negative weights are prevented.
}
\value{
  A numeric matrix where each row represents a different weight scenario. The first row contains the initial weights, 
  while subsequent rows show random modifications.
}
\seealso{
  \code{\link{gengradwei}}
}
\author{
  Cagatay Cebeci
}

\examples{
# Define initial weights
init_weights <- c(C1 = 0.4, C2 = 0.3, C3 = 0.3)

# Generate random weight modifications
scenarios <- genrandwei(init_weights, ss = 0.05, niters = 10)

# View resulting scenarios
print(scenarios)
}
