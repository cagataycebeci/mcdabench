# ORETES
oretes <- function(dmatrix, bcvec, weights, tiesmethod = "average", domplot=FALSE) {
  if (!is.matrix(dmatrix) && !is.data.frame(dmatrix)) {
    stop("Decision matrix (dmatrix) must be a matrix or a data frame.")
  }
  if (!all(bcvec %in% c(1, -1)) || length(bcvec) != ncol(dmatrix)) {
    stop("The benefit/cost vector (bcvec) must contain 1 (benefit) or -1 (cost) values.")
  }
  if (!is.numeric(weights) || length(weights) != ncol(dmatrix)) {
    stop("Weights must be a numeric vector with the same length as the number of criteria.")
  }

  n <- nrow(dmatrix)
  m <- ncol(dmatrix)

  if (is.null(rownames(dmatrix))) {
    rownames(dmatrix) <- paste0("A", 1:n)
  }
  anames <- rownames(dmatrix)

  lambdamatrix <- matrix(NA, nrow = n, ncol = m)
  for (j in 1:m) {
    minval <- min(dmatrix[, j])
    maxval <- max(dmatrix[, j])
    lambdamatrix[, j] <- if (minval == maxval) 0.5 else if (bcvec[j] == 1) (dmatrix[, j] - minval) / (maxval - minval) else (maxval - dmatrix[, j]) / (maxval - minval)
  }

  psimatrix <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) {
    psimatrix[i, ] <- rank(-lambdamatrix[i, ], ties.method = tiesmethod) 
  }

  dominance_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(anames, anames))
  for (a in 1:n) {
    for (b in 1:n) {
      if (a == b) next
      sum_weighted_preference <- sum(weights[which(lambdamatrix[a, ] > lambdamatrix[b, ])])
      dominance_matrix[a, b] <- sum_weighted_preference
    }
  }

  maxposdom <- sum(weights)
  nordommat <- if (maxposdom > 0) dominance_matrix / maxposdom else matrix(0, nrow = n, ncol = n)

  outranking_flow <- numeric(n)
  for (i in 1:n) {
    if (i > n || i <= 0) stop("Index out of bounds in dominance matrix")
    outgoing_flow <- sum(nordommat[i, ])
    incoming_flow <- sum(nordommat[, i])
    outranking_flow[i] <- outgoing_flow - incoming_flow
  }
  names(outranking_flow) <- anames

  ranking <- rank(-outranking_flow, ties.method = tiesmethod)

  results <- list(
    lambdamatrix = lambdamatrix,
    psimatrix = psimatrix,
    dominance_matrix = nordommat,
    outranking_flow = outranking_flow,
    rank = ranking
  )

  if (domplot && any(nordommat > 0)) { 
    graph_edges <- which(nordommat > 0, arr.ind = TRUE)
    edge_list <- data.frame(
      from = rownames(nordommat)[graph_edges[, 1]],
      to = colnames(nordommat)[graph_edges[, 2]],
      weight = nordommat[graph_edges]
    )
  
    g <- igraph::graph_from_data_frame(edge_list, directed = TRUE)
    igraph::plot.igraph(g, edge.label = round(igraph::E(g)$weight, 2), 
       main = "ORESTE Dominance Graph", 
       vertex.color = "lightblue", vertex.size = 30,
       edge.width = igraph::E(g)$weight * 5
    )
  } else if (domplot) {
    message("No dominance relationships found. Graph will not be plotted.")
  }
  return(results)
}