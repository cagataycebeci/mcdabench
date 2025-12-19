# Random Weight Modifications for MCDA Benchmarking
# Randomly modifies the initial weights for n times with a step size (ss)
#
genrandwei <- function(initwei, ss=NULL, niters=5) {
  m <- length(initwei)
  
  if(is.null(ss)) ss <- min(initwei) * 0.05
  
  if(ss <= 0 || ss > max(initwei)){
    stop("Step size (ss) must not be greater than maximum weight or less than 0.")
  }

  cnames <- names(initwei)
  if (is.null(cnames)) cnames <- paste0("C_", 1:m)
  
  scenario_names <- c("RW_0", paste0("RW_", 1:niters))

  scenarios <- matrix(NA, nrow=length(scenario_names), ncol=m, dimnames=list(scenario_names, cnames))

  scenarios["RW_0", ] <- initwei

  # Modify weights incrementally
  for (i in 1:niters) {
    repeat {
      modwei <- initwei + runif(length(initwei), -abs(ss), abs(ss))
      modwei <- modwei / sum(modwei)  # Normalize
      
      # Ensure no negative weights exist
      if (all(modwei > 0)) break 
    }

    # Adjust sum to exactly 1
    diff <- 1 - sum(modwei)
    modwei[which.max(modwei)] <- modwei[which.max(modwei)] + diff
    
    scenarios[i + 1, ] <- modwei
  }

  return(scenarios)
}
