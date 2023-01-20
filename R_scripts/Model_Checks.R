# Purpose: To validate triple separable GMRF models for EBS pollock 
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date 1/19/23


#' Title Returns marginal AIC for TMB models
#'
#' @param optim_model optimized model
#'
#' @return Integer of AIC values
#' @export
#'
#' @examples
#' 
margAIC <- function(optim_model) {
  
  # Get number of parameters 
  k <- length(optim_model[["par"]])
  
  # Extract objective function
  nLL <- optim_model[["objective"]]
  
  # Calculate AIC
  margAIC <- 2*k + 2*nLL
  
  return(margAIC)
}

# Set Up ------------------------------------------------------------------

# load in packages
library(here)
library(tidyverse)

# Load in models
load(here("output", "ebs_pollock_waa_models.RData")) # obj is called models

# Generate full factorial design for naming models
map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>% 
  data.frame()

# Model Checking ----------------------------------------------------------

AIC_models <- vector() # store aic
nlminb_conv <- vector() # store nlminb convergence diagnostic
grad_conv <- vector() # store maximum gradient
grad_name_conv <- vector() # store parameter with maximum gradient
hess_conv <- vector() # whether or not Hessian is positive definite 
par_values <- data.frame(rho_y = NA, rho_a = NA, 
                         rho_c = NA, log_sigma2 = NA) # parameter values
par_sd_values <- data.frame(rho_y_sd = NA, rho_a_sd = NA, 
                            rho_c_sd = NA,  log_sigma2_sd = NA) # se values

# Extract model values
for(i in 1:length(models)) {
  
  # Convergence diagnostics
  AIC_models[i] <- margAIC(models[[i]]$optim) # get aic
  nlminb_conv[i] <- models[[i]]$optim$convergence # get nlminb convergence
  grad_conv[i] <- max(abs(models[[i]]$sd_rep$gradient.fixed)) # max gradient
  # Get max gradient parameter name
  grad_name_conv[i] <- names(models[[i]]$sd_rep$par.fixed)[which.max(abs(models[[i]]$sd_rep$gradient.fixed))]
  hess_conv[i] <- models[[i]]$sd_rep$pdHess # whether or not this is pd HESS
  
  # Grab parameter values here and store
  par_values[i,] <- models[[i]]$sd_rep$par.fixed[match(colnames(par_values),
                                                       names(models[[i]]$sd_rep$par.fixed))]
  
  # Get parameter std values here and store
  par_sd_values[i,] <- sqrt(
    diag(models[[i]]$sd_rep$cov.fixed)[match(colnames(par_values),
                                             names(diag(models[[i]]$sd_rep$cov.fixed)))]
  )
  
  if(i == length(models)) {
    model_diag <- data.frame(AIC = AIC_models, nlminb_conv = nlminb_conv,
                             max_grad_name = grad_name_conv,
                             max_grad = grad_conv, pd_Hess = hess_conv)
    
    # Bind parameter values and sd together
    model_diag <- cbind(model_diag, par_values, par_sd_values)
  } # when done w/ evaluating all models
  
} # end i loop


# Munge into plot format --------------------------------------------------

# Create model names to differentiate models
model_names <- map_factorial %>% 
  mutate(rho_y_lab = case_when(rho_y == 1 ~ "y"),
         rho_a_lab = case_when(rho_a == 1 ~ "a"),
         rho_c_lab = case_when(rho_c == 1 ~ "c")) %>% 
  dplyr::select(rho_y_lab, rho_a_lab, rho_c_lab) %>% 
  dplyr::rowwise() %>% 
  tidyr::unite('model', na.rm = TRUE)

# Input model names above into model_diag df
model_diag$model <- model_names$model
# No correlation parameters estimated
model_diag$model[model_diag$model == ""] <- "None"

# Pivot this dataframe longer 
# Parmeter MLEs here only (doing ses and binding in 2 steps)
model_pars_long <- model_diag %>% 
  dplyr::select(-rho_a_sd, -rho_y_sd, -rho_c_sd, -log_sigma2_sd) %>% 
  pivot_longer(cols = c(rho_a, rho_c, rho_y, log_sigma2), 
               names_to = "parameters",  values_to = "mle_val")

# Get SE values now
model_se_long <- model_diag %>% 
  dplyr::select(rho_a_sd, rho_y_sd, rho_c_sd, log_sigma2_sd) %>% 
  pivot_longer(cols = everything(),names_to = "sd",  values_to = "sd_val")

# Now bind, these two together
model_diag_long <- cbind(model_pars_long, model_se_long)

# Create lwr and upr confidence intervals here
model_diag_long <- model_diag_long %>% 
  mutate(lwr_95 = mle_val - (1.96 * sd_val),
         upr_95 = mle_val + (1.96 * sd_val))

# Output to csv
write.csv(model_diag_long, here("output", "model_diag_vals.csv"))
