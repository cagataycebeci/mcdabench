# Normalization Functions
calcnormal <- function(X, bcvec, type = "maxmin") {
 
  # No normalization
  nullnorm <- function(x, bc = -1) {
    return(x) 
  }

  # Linear Ratio Normalization (aka Max normalization)
  rationorm <- function(x, bc = -1) { 
    if (bc == 1) {
      nvec <- x / max(x) 
    } else { 
      #nvec <- min(x) / x + 1e-6
      nvec <-  1 - x / max(x)
    } 
    return(nvec) 
  }
 
  # Max-Min Normalization
  maxminnorm <- function(x, bc = -1) {
    if (min(x) == max(x)){
      nvec <- rep(1, length(x)) # nvec <- rep(0.5, length(x))
    } else if (bc == 1) {
      nvec <- (x - min(x)) / (max(x) - min(x))
    } else {
      nvec <- (max(x) - x) / (max(x) - min(x))
    }
    return(nvec)
  }
  
  # Vector Normalization
  vecnorm <- function(x, bc = -1) {
    if (bc == 1) {
      nvec <- x / sqrt(sum(x^2) + 1e-6)
    } else {
      # nvec <- (sqrt(sum(x^2) + 1e-6)) / (x + 1e-6)
      nvec <- 1 - (x / (sqrt(sum(x^2) + 1e-6)))   
    }
    return(nvec) 
  }

  # Sum Normalization
  sumnorm <- function(x, bc = -1) {
    if (bc == 1) {
      nvec <- x / sum(x)
    } else {
      nvec <- (1 / x) / sum(1 / x)
    }
    return(nvec)
  }
   
  # Zavadskas-Turskis Normalization
  #ztnorm <- function(x, bc = -1) {
  #  if (bc == 1) {
  #    nvec <- 1 - abs((max(x) - x) / max(x))
  #  } else {
  #    nvec <- 1 - abs((min(x) - x) / min(x))
  #  }
  #  return(nvec)
  #}

  ztnorm <- function(x, bc = -1) {
    if (bc == 1) {
      nvec <- x / max(x)
    } else {
      nvec <- min(x) / x
    }
    return(nvec)
  }

  # Enhanced Accuracy Normalization
  enaccnorm <- function(x, bc = -1) {
    if (bc == 1) {
      nvec <- 1 - ((max(x) - x) / sum(max(x) - x))
    } else {
      nvec <- 1 - ((x - min(x)) / sum(x - min(x)))
    }
    return(nvec)
  }
  
  # Linear Normalization (Same with ratio)
  linnorm <- function(x, bc = -1) {
    if (!is.numeric(x)) {
      stop("Input must be a numeric vector.")
    }
    if (bc == 1) {
      nvec <- x / max(x)
    } else {
      nvec <- min(x) / (x + 1e-6)
    }
    return(nvec)
  }
  
  # Non-linear Normalization
  nonlinnorm <- function(x, bc = -1) {
    if (bc == 1) {
      nvec <- (x / max(x)) ^ 2
    } else {
      nvec <- (min(x) / x) ^ 3
    }
    return(nvec)
  }
  
  # Logaritmic Normalization
  lognorm <- function(x, bc = -1) {
    prod <- prod(x)
    if (bc == 1) {
       nvec <- log(x) / log(prod)
   } else {
      nvec <- (1 - (log(x) / log(prod))) / (length(x) - 1)
    }
    return(nvec)
  }
 
  zscore <- function(column, bc) {
    return((column - mean(column)) / sd(column))
  }

#-------------------------------------------------------------------  
  
  X <- as.matrix(X)

  colnorm <- function(column, bc, type) {
    if (type == "ratio") {
      return(rationorm(column, bc))
    } else if (type == "maxmin") {
      return(maxminnorm(column, bc))
    } else if (type == "vector") {
      return(vecnorm(column, bc))
    } else if (type == "sum") {
      return(sumnorm(column, bc))
    } else if (type == "zavadskas") {
      return(ztnorm(column, bc))
    } else if (type == "enhanced") {
      return(enaccnorm(column, bc))
    } else if (type == "linear") {
      return(linnorm(column, bc))
    } else if (type == "nonlinear") {
      return(nonlinnorm(column, bc))
    } else if (type == "logarithmic") {
      return(lognorm(column, bc))
    } else if (type == "enacc") {
      return(enaccnorm(column, bc))
    } else if (type == "zscore") {
      return(zscore(column, bc))
    } else {
      return(nullnorm(column, bc))
    }
  }
  
  if(is.null(type)) type <- "null"
  if (!is.matrix(X)) stop("X must be a valid decision matrix with dimension n x m")
  if (!all(bcvec %in% c(-1, 1))) stop("The bcvec should contain -1 for cost and 1 for benefit")
  
  normmethods <- c("enacc", "linear", "nonlinear", "logarithmic", "enhanced", "zavadskas", "sum", "vector",  "maxmin", "null", "ratio", "zscore") 
  if (!type %in% normmethods) {
    stop(paste("Invalid normalization method! The type should be set to one of: ", paste(normmethods, collapse = ", ")))
  }

  nmatrix <- matrix(NA, nrow=nrow(X), ncol=ncol(X))
  colnames(nmatrix) <- colnames(X)
  rownames(nmatrix) <- rownames(X)
  for(j in 1:length(bcvec)){
    normvec <- colnorm(X[,j], bcvec[j], type=type)
    nmatrix[,j] <- normvec
  }
  return(nmatrix)
}

