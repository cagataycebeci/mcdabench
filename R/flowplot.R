# Dominance flow plot
flowplot <- function(rankvec, colpal, txtcol = "black", orientation = "horizontal") {
  rankorder <- order(rankvec)
  if (is.null(names(rankvec))) {
    names(rankvec) <- paste0("A", seq_along(rankvec))
  }

  sorted_rankvec <- rankvec[rankorder]
  sorted_names <- names(rankvec)[rankorder]

  edges <- character()
  edgelabs <- character()
  uniqueranks <- unique(sorted_rankvec)

  rank_to_color <- setNames(colpal[seq_along(uniqueranks)], as.character(uniqueranks))
  node_colors_map <- sapply(sorted_rankvec, function(x) rank_to_color[as.character(x)])
  names(node_colors_map) <- sorted_names

  for (i in seq_len(length(sorted_rankvec) - 1)) {
    edges <- c(edges, sorted_names[i], sorted_names[i + 1])
    edgelabs <- c(edgelabs, if (sorted_rankvec[i] == sorted_rankvec[i + 1]) "=" else ">")
  }

  g <- igraph::make_empty_graph(n = length(sorted_names), directed = TRUE)
  igraph::V(g)$name <- sorted_names

  if (length(edges) > 0) {
    g <- igraph::add_edges(g, edges)
  }

  n <- length(sorted_rankvec)
  vertex_base_size <- 16
  vertex_size <- max(6, vertex_base_size - (n * 0.5))
  spacing_multiplier <- 30

  layout_func <- if (orientation == "horizontal") {
    cbind(seq_len(n) * (vertex_size * spacing_multiplier), rep(0, n))
  } else {
    cbind(rep(0, n), seq(n, 1) * (vertex_size * spacing_multiplier))
  }

  vertex_colors_for_plot <- node_colors_map[igraph::V(g)$name]

  par(mar = c(2, 2, 2, 2)) 
  graphics::plot(
    g,
    layout = layout_func,
    vertex.shape = "circle",
    vertex.size = vertex_size,
    vertex.color = vertex_colors_for_plot,
    vertex.label = igraph::V(g)$name,
    vertex.label.cex = 0.9,
    vertex.label.color = txtcol,
    edge.arrow.size = 2,
    edge.width = 2,
    edge.label = edgelabs
  )

  coords <- layout_func
  el <- igraph::as_edgelist(g)
  offset <- 5 

  for (i in seq_len(nrow(el))) {
    from <- el[i, 1]
    to <- el[i, 2]
    from_coord <- coords[which(igraph::V(g)$name == from), ]
    to_coord   <- coords[which(igraph::V(g)$name == to), ]

    mid_coord <- (from_coord + to_coord) / 2
    vec <- to_coord - from_coord
    unit_vec <- vec / sqrt(sum(vec^2))
    ort_vec <- c(-unit_vec[2], unit_vec[1])
    label_coord <- mid_coord + offset * ort_vec

    graphics::text(
      x = label_coord[1],
      y = label_coord[2],
      labels = edgelabs[i],
      col = "black", cex = 1.0, font = 2
    )
  }
}
