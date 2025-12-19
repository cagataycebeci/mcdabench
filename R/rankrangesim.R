  rankrangesim <- function (mat, k = 3, plotsim = FALSE, low_color = "blue", high_color = "red", text_color = "white") {
  
  # Function to plot heatmap
  plot_heatmap <- function(percentage_matrix, low_color = "blue", high_color = "red", text_color = "white") {
    heatmap_data <- as.data.frame(as.table(percentage_matrix))
    colnames(heatmap_data) <- c("Method1", "Method2", "Similarity_Percentage")

    p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data[["Method1"]], y = .data[["Method2"]], fill = .data[["Similarity_Percentage"]])) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = low_color, high = high_color) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "K-Range Similarity Heatmap (%)", x = "Method", y = "Method") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["Similarity_Percentage"]], 1)), color = text_color, size = 4)

    return(p)
  }
  
  # Check the inputs
  if (length(k) == 1) {
    if (k > 0) {
      k <- c(1, k)
    } else {
      k <- c(ncol(mat) + k, ncol(mat)) 
    }
  }

  if (length(k) != 2) stop("Range k must contain one or two integers to define the starting and end of range.")
  if (k[1] > k[2]) stop("Error: k1 cannot be greater than k2!")  
  if (k[1] < 1 || k[2] > ncol(mat)) stop(paste("Error: k1 must be at least 1, k2 must be at most", ncol(mat), "!"))
  if (any(k <= 0)) stop("Error: k values cannot be negative or zero!")

  results <- list()
  method_names <- rownames(mat)

  for (method in method_names) {
    row_values <- as.numeric(mat[method, ])
    sorted_indices <- order(row_values)  
    selected_alternatives <- colnames(mat)[sorted_indices[k[1]:k[2]]]  

    results[[method]] <- selected_alternatives
  }

  # Create similarity matrix
  similarity_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
  colnames(similarity_matrix) <- method_names
  rownames(similarity_matrix) <- method_names

  for (i in 1:length(method_names)) {
    for (j in 1:length(method_names)) {
      similarity_matrix[i, j] <- length(intersect(results[[method_names[i]]], results[[method_names[j]]]))
    }
  }

  percentage_matrix <- round((similarity_matrix / (k[2] - k[1] + 1)) * 100, 2)

  if (plotsim) {
    print(plot_heatmap(percentage_matrix, low_color = low_color, high_color = high_color, text_color = text_color))
  }

  return(list(selected_k = results, similarity = similarity_matrix, percentage = percentage_matrix))
}
