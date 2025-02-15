library(here)
library(tidyverse)
library(RTMB)

# Model -------------------------------------------------------------------

#' Title Constructor algorithin for correlations within ages, years, and cohort
#'
#' @param n_ages Number of ages
#' @param n_yrs Number of years
#' @param pcorr_age correlations for age
#' @param pcorr_year correaltions for year
#' @param pcorr_cohort correlaitons for cohort
#' @param ln_var_value log space variance
#' @param Var_Param variance type == 0, marginal (stationary and slower run time), == 1 conditional (non-statationary, faster run time)
#'
#' @returns Sparse precision matrix dimensioned by n_ages * n_years, n_ages * n_years
#' @export
#'
#' @examples
#' 
Get_3d_precision <- function(n_ages, n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Param){
  
  require(Matrix)
  
  index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
  i = j = x = numeric(0) # initialize posiiton to fill in precision matrix
  var_value = exp(ln_var_value) # transform to normal space
  
  for(n in 1:nrow(index)){
    age = index[n,1] # get age index out of all index combinations
    year = index[n,2] # get year index out of all index combinations
    if(age > 1 ){
      i = c(i, n)
      j = c(j, which(index[,1] == (age-1) & index[,2] == year))
      x = c(x, pcorr_year) # year correaltion indexing
    }
    if(year > 1){
      i = c(i, n)
      j = c(j, which(index[,1]==age & index[,2]==(year-1)) )
      x = c(x, pcorr_age) # age correlation indexing
    }
    if( age>1 & year>1 ){
      i = c(i, n)
      j = c(j, which(index[,1]==(age-1) & index[,2] == (year-1)) )
      x = c(x, pcorr_cohort) # cohort correlation indexing
    }
  } # end n loop
  
  # create B path matrix
  B = matrix(0, nrow = n_ages * n_yrs, ncol = n_ages * n_yrs)
  B[cbind(i, j)] = x
  B = as(B, "sparseMatrix") 
  
  # identity matrix
  I = as(diag(1, n_ages * n_yrs, n_ages * n_yrs), "sparseMatrix")
  
  if(Var_Param == 0) d = var_value # conditional variance (non-stationary variance)
  
  # Solve Omega recursively for stationary variance (accumulator function)
  if(Var_Param == 1) {
    L = solve(I-B) # solve to get accumulator function for stationary variance
    d = rep(0, nrow(index))
    for(n in 1:nrow(index) ){
      if(n==1){
        d[n] = var_value
      }else{
        cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
        d[n] = (var_value-cumvar) / L[n,n]^2
      }
    } # end n loop
  } # end marginal variance (stationary variance)
  
  # omega matrix
  Omega_inv = diag(1/d, n_ages * n_yrs, n_ages * n_yrs)
  Q = as((I-t(B)) %*% Omega_inv %*% (I-B), "sparseMatrix") # solve for precision
  
  return(Q)
}

growth_3d = function(pars) {
  
  RTMB::getAll(pars, data) # load in starting values and data

  # make precision
  Q_sparse = Get_3d_precision(nrow(ln_Y_at), ncol(ln_Y_at), rho_a, rho_y, rho_c, log_sigma2, Var_Param)
  
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
  jnLL = jnLL - dgmrf(x = as.vector(eps_at), mu = 0, Q = Q_sparse, TRUE)
  
  RTMB::REPORT(jnLL);
  RTMB::REPORT(Q_sparse);
  RTMB::REPORT(mu_at);
  RTMB::REPORT(ln_Y_at); 
  
  return(jnLL)
}


# Run Model ---------------------------------------------------------------

# Load in WAA matrix (only use fishery data)
waa_df <- read.csv(here("data", "ebs_waa.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Load in std for WAA matrix
waa_std_df <- read.csv(here("data", "ebs_waa_std.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Number of projection years
n_proj_years <- 5

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

ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), 
                                  "year" = seq_len(length(years) + n_proj_years) ))

# Now, input these components into a data list
data <- list( years = years,
              ages = ages,
              X_at = X_at,
              Xsd_at = Xsd_at,
              ay_Index = ay_Index,
              n_proj_years = n_proj_years,
              Var_Param = 0) # conditional

# Input parameters into a list
parameters <- list( rho_y = 0,
                    rho_a = 0,
                    rho_c = 0,
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

start_time = Sys.time()
# make AD model function
growth_3d_model <- RTMB::MakeADFun(growth_3d, parameters = parameters, map = map, random = 'ln_Y_at')

# Now, optimize the function
growth_optim <- stats::nlminb(growth_3d_model$par, growth_3d_model$fn, growth_3d_model$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

# newton steps
try_improve <- tryCatch(expr =
                          for(i in 1:3) {
                            g = as.numeric(growth_3d_model$gr(growth_optim$par))
                            h = optimHess(growth_optim$par, fn = growth_3d_model$fn, gr = growth_3d_model$gr)
                            growth_optim$par = growth_optim$par - solve(h,g)
                            growth_optim$objective = growth_3d_model$fn(growth_optim$par)
                          }
                        , error = function(e){e}, warning = function(w){w})

max(growth_3d_model$gr())
growth_3d_model$optim <- growth_optim # Save optimized model results
growth_3d_model$sd_rep <- RTMB::sdreport(growth_3d_model) # Get sd report
end_time = Sys.time()
end_time - start_time
growth_3d_model$rep <- growth_3d_model$report(growth_3d_model$env$last.par.best) # Get report
growth_3d_model$sd_rep

reshape2::melt(exp(growth_3d_model$rep$ln_Y_at)) %>% 
  ggplot(aes(x = Var2, y = value)) +
  geom_line() +
  facet_wrap(~Var1, scales = 'free')


