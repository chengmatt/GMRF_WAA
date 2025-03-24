growth_2d = function(pars) {
  
  RTMB::getAll(pars, data) # load in starting values and data
  
  # Function to constrain values between -1 and 1
  rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1
  
  jnLL = 0 # set up jnLL 
  
  # set up containers
  mu_at = array(0, dim = c(nrow(ln_Y_at), ncol(ln_Y_at)))
  eps_at = array(0, dim = c(nrow(ln_Y_at), ncol(ln_Y_at)))
  
  # transform pars
  L0 = exp(ln_L0);
  Linf = exp(ln_Linf);
  k = exp(ln_k);
  alpha = exp(ln_alpha);
  beta = exp(ln_beta);
  
  # get parametric form
  for(a in 1:nrow(X_at)) {
    for(t in 1:ncol(X_at)) {
      mu_at[a,t] = Linf - (Linf - L0) * exp(-k * a)
      mu_at[a,t] = alpha * mu_at[a,t]^beta
    }
  }
  
  # observation likelihood
  for(a in 1:nrow(X_at)) {
    for(t in 1:ncol(X_at)) {
      if(!is.na(X_at[a,t])) {
        if(Xsd_at[a,t] > 0) jnLL = jnLL + -dnorm(log(X_at[a,t]), ln_Y_at[a,t], Xsd_at[a,t], TRUE)
      }
    }
  }
  
  # process error liklelihood
  eps_at = ln_Y_at - log(mu_at)
  rho_a_trans = rho_trans(rho_a) # correaltion across ages
  rho_y_trans = rho_trans(rho_y) # correaltion across ages
  sigma2 = exp(log_sigma2)^2 # get sigma
  
  if(type == 0) { # 2d ar1
    unit_var = sqrt(sigma2) / sqrt(1 - rho_y_trans^2) / sqrt(1 - rho_a_trans^2) # Define 2d unit variance
    # Define ar1 separable functions
    f1 = function(x) dautoreg(x, mu = 0, phi = rho_a_trans, log = TRUE)
    f2 = function(x) dautoreg(x, mu = 0, phi = rho_y_trans, log = TRUE)
    jnLL = jnLL - RTMB::dseparable(f1, f2)(eps_at, scale = unit_var)
  }
  
  if(type == 1) { # iid along ages, 1dar1 along years
    f1 = function(x) sum(dnorm(x, 0, sqrt(sigma2), log = TRUE))
    f2 = function(x) dautoreg(x, mu = 0, phi = rho_y_trans, log = TRUE)
    jnLL = jnLL - RTMB::dseparable(f1, f2)(eps_at)
  }
  
  # Quick notes about dseparable:
  # Can take any number of log densities and treats them as separable to formulate a joint MVN
  # Separable log densities must be mean zero.
  # dseparable outputs log density by default.
  # Expects either a matrix or array input for x; can't be a vector
  # In this case, f1 acts on the age row, f2 acts on the time column.
  # If eps_at was an array, f1 would act on the first slice, f2 would be on the second slice, and f3 would be on the third slice
  
  # So if we have type == 1, and we switched the ordering of f2 and f1, it would apply ar1 on ages, but iid on years (i.e., ordering matters)
  
  RTMB::REPORT(jnLL);
  RTMB::REPORT(mu_at);
  RTMB::REPORT(ln_Y_at); 
  
  return(jnLL)
}

# 2d model ----------------------------------------------------------------
# Load in WAA matrix (only use fishery data)
waa_df <- read.csv(here("data", "ebs_waa.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Load in std for WAA matrix
waa_std_df <- read.csv(here("data", "ebs_waa_std.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Number of projection years
n_proj_years <- 30

# Years
years <- waa_df$year

# Ages (goes from age 3 - 15+)
ages <- parse_number(colnames(waa_df)[-1])

# Read in data weight at age matrix
X_at <- t(as.matrix(waa_df[,-1])) # removing first col (year column)

# Create projection columns
proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years)

# Append NA for projection year
X_at <- cbind(X_at, proj_cols) 

# Read in standard deviations for weight at age matrix
Xse_at <- t(as.matrix(waa_std_df[,c(-1)])) # removing first col (year column) 

# Convert to CV
Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )

# Now convert back to sd in lognormal space
Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))

# Now, input these components into a data list
data <- list( years = years,
              ages = ages,
              X_at = X_at,
              Xsd_at = Xsd_at,
              n_proj_years = n_proj_years,
              type = 1) 

# Input parameters into a list
parameters <- list( rho_y = 0,
                    rho_a = 0,
                    log_sigma2 = log(0.1),
                    ln_L0 = log(45),
                    ln_Linf = log(80),  # Fixed at arbitrary value
                    ln_k = log(0.15),
                    ln_alpha = log(3.5e-7), # Start alpha at a reasonable space 
                    # Starting value for alpha derived from a run where none of the rhos were estimated.
                    ln_beta = log(3), # Fix at isometric
                    ln_Y_at = array(0,dim=dim(X_at)) ) 

map <- list(
  ln_Linf = factor(NA), ln_beta = factor(NA)
)

# turn rho_a off if age is assumed to be iid
if(data$type == 1) map$rho_a = factor(NA)

start_time = Sys.time()
# make AD model function
growth_2d_model <- RTMB::MakeADFun(growth_2d, parameters = parameters, map = map, random = 'ln_Y_at')

# Now, optimize the function
growth_optim <- stats::nlminb(growth_2d_model$par, growth_2d_model$fn, growth_2d_model$gr,
                              control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

# newton steps
try_improve <- tryCatch(expr =
                          for(i in 1:3) {
                            g = as.numeric(growth_2d_model$gr(growth_optim$par))
                            h = optimHess(growth_optim$par, fn = growth_2d_model$fn, gr = growth_2d_model$gr)
                            growth_optim$par = growth_optim$par - solve(h,g)
                            growth_optim$objective = growth_2d_model$fn(growth_optim$par)
                          }
                        , error = function(e){e}, warning = function(w){w})

max(growth_2d_model$gr())
growth_2d_model$optim <- growth_optim # Save optimized model results
growth_2d_model$sd_rep <- RTMB::sdreport(growth_2d_model) # Get sd report
end_time = Sys.time()
end_time - start_time
growth_2d_model$rep <- growth_2d_model$report(growth_2d_model$env$last.par.best) # Get report
growth_2d_model$sd_rep



# Model Comparison --------------------------------------------------------

reshape2::melt(exp(growth_2d_model$rep$ln_Y_at)) %>% mutate(type = '2d') %>% 
  mutate(cohort = Var2 - Var1) %>% 
  ggplot(aes(x = Var2, y = value, color = type)) +
  geom_line() +
  facet_wrap(~Var1, scales = 'free')
