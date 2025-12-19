# Weights Calculation for Multi-criteria Decision Making
#
calcweights <- function(dmatrix, bcvec=NULL, type="equal", normethod=NULL){

    X <- as.matrix(dmatrix)  

    # Equal weight method
    equalwei <- function(mat) {
       n <- ncol(mat)
       return(rep(1 / n, n))
    }

    # Standard deviation weight method
    sdwei <- function(mat) {
      sdev <- apply(mat, 2, sd)
      return(sdev / sum(sdev))      
    }

    # Variance weight method
    varwei <- function(mat) {
      vars <- apply(mat, 2, var)
      return(vars / sum(vars))      
    }

    # Critic weighting method
    critwei <- function(mat) {     
      sdev <- apply(mat, 2, sd)
      coef <- cor(mat, method = "pearson")
      C <- sdev * colSums(1 - coef)
      return(C / sum(C))
    }
    
    # Entropy method
    entwei <- function(mat) {
      m <- nrow(mat)
      n <- ncol(mat)
    
      col_sums <- colSums(mat)
      p_matrix <- mat
      for (j in 1:n) {
        if (col_sums[j] == 0) {
          p_matrix[, j] <- 0
        } else {
          p_matrix[, j] <- mat[, j] / col_sums[j]
        }
      }
      
      entropy_terms <- p_matrix * log(p_matrix)
      entropy_terms[is.nan(entropy_terms)] <- 0
      entropy_terms[is.infinite(entropy_terms)] <- 0
    
      sum_p_ln_p <- colSums(entropy_terms)
      
      ln_m <- log(m)
      
      if (ln_m == 0) {
        entropies <- rep(0, n)
      } else {
        entropies <- -sum_p_ln_p / ln_m
      }
    
      E_diffs <- 1 - entropies
      sum_of_E_diffs <- sum(E_diffs)
    
      if (sum_of_E_diffs == 0) {
        warning("Weights are equally distributed because all entropy differences (1 - E_j) are zero.")
        weights <- rep(1/n, n)
      } else {
        weights <- E_diffs / sum_of_E_diffs
      }
    
      if (!is.null(colnames(mat))) {
        names(weights) <- colnames(mat)
      }
    
      return(weights)
    }
    
    entropy <- function(x_matrix) {
      m <- nrow(x_matrix)
      n <- ncol(x_matrix)
    
      p_matrix <- matrix(NA, nrow = m, ncol = n)
      for (j in 1:n) {
        col_sum <- sum(x_matrix[, j])
        p_matrix[, j] <- x_matrix[, j] / col_sum
      }
    
      E_j_vector <- numeric(n)
      for (j in 1:n) {
        term_sum <- 0
        for (i in 1:m) {
          if (p_matrix[i, j] > 0) {
            term_sum <- term_sum + (p_matrix[i, j] * log(p_matrix[i, j]))
          }
        }
        E_j_vector[j] <- -term_sum / log(m)
      }
    
      w_j_vector <- numeric(n)
      sum_of_1_minus_E_j <- sum(1 - E_j_vector)
      for (j in 1:n) {
        w_j_vector[j] <- (1 - E_j_vector[j]) / sum_of_1_minus_E_j
      }
    
      names(w_j_vector) <- colnames(x_matrix)
      return(w_j_vector)
}


    # Alternative Entropy method
    entwei2 <- function(mat) {
      m <- nrow(mat)
      n <- ncol(mat)
      pij <- matrix(0, nrow = m, ncol = n)
    
      for (j in 1:n) {
        if (sum(mat[, j]) != 0) {
          pij[, j] <- mat[, j] / sum(mat[, j])
        }
      }
    
      ej <- -colSums(pij * log(pij))
      ej[is.nan(ej)] <- 0  # Replace NaNs with 0
      weights <- (1 - ej) / sum(1 - ej)
      return(weights)
    }

    # GINI weighting method
    giniwei <- function(mat) {     
      n <- nrow(mat)
      m <- ncol(mat)
      weights <- numeric(m)
      
      for (j in 1:m) {
        values <- numeric(n)
        for (i in 1:n) {
          values[i] <- sum(abs(mat[i, j] - mat[, j]) / 
           (2 * n^2 * (sum(mat[, j]) / n)))
        }
        weights[j] <- sum(values)
      }
      
      return(weights / sum(weights))
    }
    
    # MEREC weighting method
    merecwei <- function(dmatrix, bcvec) {
      n <- nrow(dmatrix)
      m <- ncol(dmatrix)

      if (length(bcvec) != m) {
        stop("Error: The length of bcvec must match the number of columns in decision matrix.")
      }

      maxminorm <- function(data, bc) {
        min_val <- min(data)
        max_val <- max(data)
        if (bc == 1) {
          return((data - min_val) / (max_val - min_val))
        } else if (bc == -1) {
          return((max_val - data) / (max_val - min_val))
        } else {
          stop("Error: Invalid bcvec value. The values must be -1 (cost) or 1 (benefit).")
        }
      }
      mat <- matrix(0, nrow = n, ncol = m)
      for (j in 1:m) {
        mat[, j] <- maxminorm(dmatrix[, j], bcvec[j])
      }
      epsilon <- .Machine$double.eps
      S <- log(1 + (1 / n * colSums(abs(log(mat + epsilon)))))
      Sp <- numeric(m)
      for (j in 1:m) {
        exmatrix <- mat[, -j, drop = FALSE]
        if (ncol(exmatrix) > 0) {
          tempval <- log(1 + (1 / n * sum(abs(log(exmatrix + epsilon)))))
          Sp[j] <- tempval
        } else {
          Sp[j] <- 0
        }
      }
      E <- abs(Sp - S)
      weights <- E / sum(E)
      return(weights)
    }
    
    # Geometric mean method
    geomwei <- function(xmatrix, bcvec) {
      n <- nrow(xmatrix)
      m <- ncol(xmatrix)

      if (length(bcvec) != m) {
        stop("Error: The length of benefit-cost vector must match the number of columns in decision matrix.")
      }

      maxminorm <- function(xdata, bc) {
        min_val <- min(xdata)
        max_val <- max(xdata)
        if (bc == 1) {
          return((xdata - min_val) / (max_val - min_val))
        } else if (bc == -1) {
          return((max_val - xdata) / (max_val - min_val))
        } else {
          stop("Error: Invalid value in benefic-cost vector. The values must be -1  for cost or 1  for benefit.")
        }
      }
          
      mat <- matrix(0, nrow = n, ncol = m)
      for (j in 1:m) {
        mat[, j] <- maxminorm(xmatrix[, j], bcvec[j])
      }

      geometric_means <- apply(mat + .Machine$double.eps, 2, function(x) {
        prod(x)^(1/n)
      })

      weights <- geometric_means / sum(geometric_means)
      return(weights)
    }
    
    # ROC method
    rocwei <- function(mat) {
      n <- ncol(mat)
      kvals <- 1:n
      weights <- (1 / kvals) / sum(1 / kvals)
      return(weights)
    }
    
    # RS weight
    rswei <- function(mat) {
      n <- ncol(mat)
      j <- 1:n
      weights <- 2 * (n + 1 - j) / (n * (n + 1)) 
      return(weights)
    }
    
    # Modified Preference Selection Index (MPSI)
    mpsiwei <- function(mat) {
      n <- nrow(mat)
      m <- ncol(mat)  
      nmatrix <- matrix(NA, nrow = n, ncol = m)
    
      # Maxmin normalization
      for (j in 1:m) {
        min_val <- min(mat[, j])
        max_val <- max(mat[, j])
        nmatrix[, j] <- (mat[, j] - min_val) / (max_val - min_val)
      }   
    
      cmeans <- colMeans(nmatrix)
    
      prefvar <- numeric(m)
      for (j in 1:m) { 
        prefvar[j] <- sum((nmatrix[, j] - cmeans[j])^2)
      }
    
      return(prefvar / sum(prefvar))
    }

    # Angle weighting
    anglewei <- function(mat) {
        
      if (ncol(mat) == 1) {
        return(c(1))  # If there is only one criterion, assign a weight of 1
      }
    
      nmatrix <- apply(mat, 2, function(x) {
        if (max(x) - min(x) == 0) {
          return(rep(0.5, length(x)))  # Avoid division by zero in normalization
        } else {
          return((x - min(x)) / (max(x) - min(x)))  # Min-max normalization
        }
      })
      
      nmatrix <- as.matrix(nmatrix)  # Ensure the normalized data remains a matrix
      
      # Compute cosine similarity between criteria
      simmatrix <- t(nmatrix) %*% nmatrix
      
      # Compute norms (vector lengths) for each criterion
      column_norms <- sqrt(colSums(nmatrix^2))
      
      # Display a warning if any column has a zero norm, as this may affect angle calculations
      if (any(column_norms == 0)) {
          warning("One or more criteria columns have a zero norm after normalization. This may impact angle calculations.")
      }
    
      # Compute the outer product of norms for cosine similarity denominator
      norms_matrix <- outer(column_norms, column_norms, "*")
      
      # Safely compute the angle matrix
      # Handle division by zero and restrict values to the valid range for the acos function
      valid_ratios <- ifelse(norms_matrix == 0, 0, simmatrix / norms_matrix)
      
      # Prevent floating-point precision errors by clamping values between [-1, 1]
      angmatrix <- acos(pmin(pmax(valid_ratios, -1), 1))
      
      # Replace any NaN values with zero (e.g., when calculating angles of completely zero columns)
      angmatrix[is.na(angmatrix)] <- 0
      
      # Compute raw weights by summing angle differences per criterion
      # Each weight represents how 'distinct' or 'varied' the criterion is compared to others
      raw_weights <- colSums(angmatrix)
    
      # Ensure the weights sum to 1
      if (sum(raw_weights) > 0) {
        weights <- raw_weights / sum(raw_weights)  # Normalize weights
      } else {
        # If all raw weights sum to zero (e.g., all criteria are identical or all angles are zero),
        # Assign equal weights to avoid division errors and provide a valid output
        weights <- rep(1 / length(raw_weights), length(raw_weights))
        warning("All computed raw weights sum to zero. Assigning equal weights.")
      }
      
      return(weights)
   }

#-------------------------------------------------------------------------------   
    
    # Check the arguments
    if (!is.matrix(X) || !is.numeric(X)) {
       stop("Decision matrix must be a numeric matrix.")
    }

    # Normalize if required
    if(!is.null(normethod)){
       if (is.null(bcvec)){
         stop("The bcvec, benefit-cost vector must be input if normalization is required")
       }
       X <- calcnormal(X, bcvec=bcvec, type=normethod)
    }
    
    # Weight calculation
    wmethods <- c("angle", "critic", "entropy", "entropy2", "equal", "geom", "gini", "merec", "mpsi", "roc","rs", "sdev", "var")
    if (!any(type %in% wmethods)) {
      stop(paste("Invalid method name for calculation of weights! Please input one of these:", paste(wmethods, collapse = ", ")))
    }

    result <- switch(type,
              equal = equalwei(X),
              angle = anglewei(X),
              critic = critwei(X),
              entropy = entwei(X),
              entropy2 = entwei2(X),
              geom = geomwei(X, bcvec=bcvec),
              gini = giniwei(X),
              merec = merecwei(X, bcvec=bcvec),
              mpsi = mpsiwei(X),
              roc = rocwei(X),
              rs = rswei(X),
              sdev = sdwei(X),
              var = varwei(X))

    names(result) <- colnames(X)
    return(result)
}  
