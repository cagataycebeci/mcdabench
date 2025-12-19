ranksrd <- function(rankmat, refopt = "mean", defrank = NULL, tiesmethod="average") {
  if (!is.matrix(rankmat) && !is.data.frame(rankmat)) {
    stop("rankmat must be a matrix or data frame.")
  }
  if (any(is.na(rankmat))) {
    stop("rankmat cannot contain NA values.")
  }
  
  valid_refopts <- c("mean", "best", "worst", "custom")
  if (!refopt %in% valid_refopts) {
    stop(paste("refopt must be one of:", paste(valid_refopts, collapse = ", ")))
  }

  n <- nrow(rankmat)
  anames <- rownames(rankmat)

  if (is.null(anames)) {
    stop("rankmat must have row names to identify alternatives.")
  }
  
  # Construct the Reference Ranking
  if (refopt == "mean") {
    refrank <- apply(rankmat, 1, mean)
  } else if (refopt == "best") {
    refrank <- rep(1, n)
    names(refrank) <- anames
  } else if (refopt == "worst") {
    max_rank_possible <- max(rankmat)
    refrank <- rep(max_rank_possible, n)
    names(refrank) <- anames
  } else if (refopt == "custom") {
    if (is.null(defrank)) {
      stop("For refopt = 'custom', defrank must be provided.")
    }
    if (!is.numeric(defrank)) {
      stop("defrank must be a numeric vector.")
    }
    if (length(defrank) != n) {
      stop(paste("defrank must have the same length as the number of alternatives (", n, ").", sep = ""))
    }
    # Ensure defrank has names matching anames for robust matching
    if (is.null(names(defrank))) {
      warning("defrank has no names. Assuming order matches rankmat row names.")
      names(defrank) <- anames
    } else {
      # Reorder defrank to match rankmat row names if needed
      defrank <- defrank[anames]
      if (any(is.na(defrank))) {
          stop("Names in defrank do not fully match rankmat row names.")
      }
    }
    refrank <- defrank
  }

  # Ensure refrank has names corresponding to alternative names
  # This is crucial for correct matching in the next step
  if (is.null(names(refrank))) {
    names(refrank) <- anames
  }

  srdvals <- numeric(n)
  names(srdvals) <- anames

  for (i in 1:n) {
     aname <- anames[i]
     srdvals[aname] <- sum(abs(rankmat[aname, ] - refrank[aname]))
  }

  # Rank alternatives by SRD values (from lowest to highest: best performance) 
  srdrank <- sort(srdvals)
  ranking <- rank(-refrank, ties.method=tiesmethod)

  return(list(
    srdvals = srdvals,
    srdrank = srdrank,
    refrank = refrank,
    rank = ranking
  ))
}
