\name{rankrangesim}
\alias{rankrangesim}
\encoding{UTF-8}
\title{Calculate Similarity of Ranked Alternatives in a Given Range}
\description{
Calculates the similarity between ranking methods by selecting alternatives from a given range and comparing them across different methods.
Additionally, it provides an option to visualize the similarity matrix as a heatmap.
}
\usage{
rankrangesim(mat, k = 3, plotsim = FALSE, 
   low_color = "blue", high_color = "red", text_color = "white")
}
\arguments{
  \item{mat}{A matrix containing ranked alternatives with methods as row names and alternatives as column names.}
  \item{k}{A numeric value or a vector of length 2 specifying the range of rankings to compare. If a single number is provided, 
  it selects from the beginning or end of the rankings.}
  \item{plotsim}{Logical. If TRUE, generates a heatmap showing the similarity percentages across methods.}
  \item{low_color}{Color for low similarity values in the heatmap. Default is "blue".}
  \item{high_color}{Color for high similarity values in the heatmap. Default is "red".}
  \item{text_color}{Color for the text labels inside the heatmap cells. Default is "white".}
}
\details{
This function takes a ranking matrix and extracts the top or bottom alternatives for each ranking method within the specified range.
It then calculates how often methods select the same alternatives, forming a similarity matrix.

If \code{plotsim = TRUE}, the function generates a heatmap using ggplot2 to visually display the similarity percentages.
}
\value{
A list containing:
  \item{selected_k}{A list of selected alternatives for each ranking method.}
  \item{similarity}{A matrix showing the count of common selected alternatives between methods.}
  \item{percentage}{A matrix representing similarity percentages of selected alternatives between methods.}
}
\seealso{
  \code{\link{rankcompare}}, \code{\link{rankentboot}},\code{\link{rankentperm}},
  \code{\link{rankmia}}, \code{\link{rankwilcox}}, \code{\link{rankwssim}}, 
  \code{\link{rankspearman}}
}

\examples{
# Example ranking matrix
rankmat <- matrix(c(
   3, 12, 1, 8, 4, 7, 6, 10, 2, 5, 11, 9,
   3, 12, 2, 8, 4, 9, 6, 10, 1, 5, 11, 7,
   3, 11.5, 4, 8, 1.5, 7, 5, 9.5, 1.5, 6, 11.5, 9.5,
   4, 3, 10, 7, 2, 6, 1, 9, 5, 11, 8, 12,
   6, 3, 10, 5, 2, 7, 1, 9, 4, 11, 8, 12,
   3, 12, 4, 8, 1, 7, 6, 9, 2, 5, 10, 11,
   3, 12, 1, 8, 4, 7, 6, 10, 2, 5, 11, 9,
   5, 3, 11, 7, 1, 6, 2, 9, 4, 10, 8, 12,
   7, 1, 11, 9, 2, 4, 3, 8, 6, 10, 5, 12,
   3, 9, 6, 8, 1, 5, 4, 10, 2, 7, 11, 12
), nrow = 10, byrow = TRUE)

rownames(rankmat) <- c("ARAS", "EDAS", "ELECT4", "FUCA", "GRA", 
   "MABAC", "MARCOS", "MEGAN", "MOORA", "PROMT2")
colnames(rankmat) <- paste0("G", 1:12)

# Similarity for the best and worst k
k1 <- 3 # Best 3 alternatives 
k2 <- -3 # Worst 3 alternatives 

res_k1 <- rankrangesim(rankmat, k=k1, plotsim = TRUE)  # Best 3 
res_k2 <- rankrangesim(rankmat, k=k2, plotsim = FALSE) # Worst 3 

print(res_k1$selected_k) 
print(res_k1$similarity) 
print(res_k1$percentage) 

print(res_k2$selected_k) 
print(res_k2$similarity) 
print(res_k2$percentage) 

# Similarity for a specific range 

kr <- c(5, 8) 
res_kr <- rankrangesim(rankmat, k=kr, plotsim = TRUE, 
   low_color = "dodgerblue", high_color = "orange", text_color = "gray")

print(res_kr$selected_k) 
print(res_kr$similarity) 
print(res_kr$percentage)
}
