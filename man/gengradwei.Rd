\name{gengradwei}
\alias{gengradwei}
\encoding{UTF-8}
\title{Gradual Modification of Initial Weights}
\description{Generates a matrix of weights by gradually modifying one criterion at a time while rescaling others accordingly.}
\usage{
gengradwei(initwei, rp = seq(0.1, 0.5, 0.1))
}
\arguments{
  \item{initwei}{A numeric vector of initial weights assigned to each criterion.}
  \item{rp}{A numeric sequence defining the percentage reduction applied to each criterion's weight, with default values ranging from 10\% to 50\%.}
}
\value{A numeric matrix where each row represents a different weight scenario. The first row contains the initial weights, while subsequent rows show modified weights based on the reduction percentages.}
\details{This function iteratively modifies the initial weights, reducing one criterion at a time while proportionally adjusting others to maintain a valid weight distribution. Negative weights are prevented by setting a lower bound. The resulting matrix provides different weighting scenarios, which can be useful for sensitivity analysis in multi-criteria decision-making models.}
\references{
Stevic, Z., Pamucar, D., Puška, A., & Chatterjee, P. (2020). Sustainable supplier selection in healthcare industries using a new MCDM method: Measurement of alternatives and ranking according to COmpromise solution (MARCOS). \emph{Computers & Industrial Engineering, 140}, 106231. \doi{10.1016/j.cie.2019.106231}

Demir, G., Chatterjee, P., & Pamucar, D. (2024). Sensitivity analysis in multi-criteria decision making: A state-of-the-art research perspective using bibliometric analysis. \emph{Expert Systems with Applications, 237}, 121660. \doi{10.1016/j.eswa.2023.121660}

Rachman, A. P., Ichwania, C., Mangkuto, R. A., Pradipta, J., Koerniawan, M. D., & Sarwono, J. (2024). Comparison of multi-criteria decision-making methods for selection of optimum passive design strategy. \emph{Energy and Buildings, 314}, 114285. \doi{10.1016/j.enbuild.2024.114285}
}
\seealso{\code{\link{genrandwei}}}
\author{Cagatay Cebeci}
\examples{
# Define initial weights
init_weights <- c(C1 = 0.5, C2 = 0.3, C3 = 0.2)

# Generate gradual weight modifications
scenarios <- gengradwei(init_weights)
print(scenarios)
}
