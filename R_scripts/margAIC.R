# Purpose: Function to compute marginal AIC for TMB models
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 1/16/23

#' Title Returns marginal AIC for TMB models
#'
#' @param optim_model optimized model
#'
#' @return Integer of AIC values
#' @export
#'
#' @examples
margAIC <- function(optim_model) {
  
  # Get number of parameters 
  k <- length(optim_model[["par"]])
  
  # Extract objective function
  nLL <- optim_model[["objective"]]
  
  # Calculate AIC
  margAIC <- 2*k + 2*nLL
  
  return(margAIC)
}

