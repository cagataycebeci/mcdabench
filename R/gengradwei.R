# Gradually modification of the initial weights
#
gengradwei <- function(initwei, rp = seq(0.1, 0.5, 0.1)) {

  m <- length(initwei)
  cnames <- names(initwei)
  if (is.null(cnames)) cnames <- paste0("C_", 1:m)
  
  # Create scenario names
  scenario_names <- c("init", paste0(rep(cnames, each=length(rp)), "-", round(rep(rp, times=m) * 100), "%"))

  # Initialize matrix (row names must match exactly!)
  scenarios <- matrix(NA, nrow = length(scenario_names), ncol = m, dimnames = list(scenario_names, cnames))

  # Initial scenario
  scenarios["init", ] <- initwei

  # Iterate through each criterion
  row_index <- 2
  for (i in 1:m) {
    for (j in seq_along(rp)) {
      curwei <- initwei
      W_alpha <- initwei[i]
      W_n_alpha <- W_alpha * (1 - rp[j])

      # Prevent negative weights
      if (W_n_alpha < 0) {
        W_n_alpha <- 1e-6
        warning(paste("Criterion", cnames[i], "reduced to 1e-16 due to excessive reduction percentage."))
      }

      # Calculate total remaining weights and rescale them
      W_beta <- sum(initwei[-i])
      scaling_factor <- (1 - W_n_alpha) / W_beta
      W_n_beta <- initwei[-i] * scaling_factor

      # Update weights
      curwei[i] <- W_n_alpha
      curwei[-i] <- W_n_beta

      scenarios[row_index, ] <- curwei
      row_index <- row_index + 1
    }
  }

  return(scenarios)
}
