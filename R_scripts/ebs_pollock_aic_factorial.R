# Purpose: To run a triple separable model testing for evidence of
# cohort/age/year effects for EBS pollock
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 1/16/23


# Set up ------------------------------------------------------------------

source(here("R_scripts", "margAIC.R"))

library(here)
library(tidyverse)
library(TMB)
library(cowplot)

# Compile and load in model
setwd(here("src"))
compile("triple_sep_waa.cpp")
dyn.load(dynlib("triple_sep_waa"))

# Load in WAA matrix (only use fishery data)
waa_df <- read.csv(here("data", "ebs_waa.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Load in std for WAA matrix
waa_std_df <- read.csv(here("data", "ebs_waa_std.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)


# Set up TMB data ----------------------------------------

# Years
years <- waa_df$year

# Ages (goes from age 3 - 15+)
ages <- parse_number(colnames(waa_df)[-1])

# Read in data weight at age matrix
WAA <- as.vector(as.matrix(waa_df[,-1])) # removing first col (year column)

# Read in standard deviations for weight at age matrix
WAA_std <- as.vector(as.matrix(waa_std_df[,-1])) # removing first col (year column)

# Create an index for ages and years to feed into TMB, which helps construct the precision matrix
ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), 
                                  "year" = seq_len(length(years)) ))


# Set up TMB Model --------------------------------------------------------

# Now, input these components into a data list
data <- list(years = years, ages = ages,
             WAA = WAA, WAA_std = WAA_std,  
             ay_Index = ay_Index)

# Input parameters into a list
parameters <- list(rho_y = 0, rho_a = 0,
                   rho_c = 0, log_sigma2 = rbeta(1, 1, 1),
                   WAA_re = rnorm(length(as.vector(WAA)), 1, 0.1))


# Run factorial models ----------------------------------------------------

# Generate full factorial design
map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>% 
  data.frame()

# Define number of extra newton steps we want to take
n.newton <- 3

# Empty list to store model objects
models <- list()

for(n_fact in 1:nrow(map_factorial)) {
  
  # Create empty map list object
  map <- list()
  
  # Extract combinations of parameters estimated here
  map_fact <- map_factorial[n_fact,]
  
  # Create our mapping list to feed into MakeADFun
  for(m in 1:length(names(map_fact))) {
    
    if(map_fact[1,names(map_fact)[m]] == 0) { # if factorial = 0, turn estimation off
      map[[m]] <- factor(NA)
      names(map)[m] <- names(map_fact)[m] # Name list object
    } else{ # factorail == 1
      map[[m]] <- factor(1)
      names(map)[m] <- names(map_fact)[m] # Name list object
    } # ifelse statement for constructing map list object
  } # end m loop

  # Now, make AD model function
  waa_model <- MakeADFun(data = data, parameters = parameters,
                         map = map, random = "WAA_re",
                         DLL = "triple_sep_waa")
  
  # Now, optimize the function
  waa_optim <- stats::nlminb(waa_model$par, waa_model$fn, waa_model$gr,  
                             control = list(iter.max = 1e5, eval.max = 1e5))
  
  # Take some additional newton steps to make sure we reach a minimum
  tryCatch(expr = for(i in 1:n.newton) {
    g = as.numeric(waa_model$gr(waa_optim$par))
    h = optimHess(waa_optim$par, fn = waa_model$fn, gr = waa_model$gr)
    waa_optim$par = waa_optim$par - solve(h,g)
    waa_optim$objective = waa_model$fn(waa_optim$par)
  }, error = function(e){e})
  
  # Save optimized model results
  waa_model$optim <- waa_optim
  
  # Get report
  waa_model$rep <- waa_model$report(waa_model$env$last.par.best)
  
  # Get sd report
  waa_model$sd_rep <- sdreport(waa_model)
  
  models[[n_fact]] <- waa_model

} # loop through to run multiple models



# Model Checking ----------------------------------------------------------

AIC_models <- vector() # store aic
nlminb_conv <- vector() # store nlminb convergence diagnostic
grad_conv <- vector() # store maximum gradient
grad_name_conv <- vector() # store parameter with maximum gradient
hess_conv <- vector() # whether or not Hessian is positive definite 
par_values <- data.frame(rho_y = NA, rho_a = NA, rho_c = NA, log_sigma2 = NA) # parameter values
par_sd_values <- data.frame(rho_y_sd = NA, rho_a_sd = NA, rho_c_sd = NA,  log_sigma2_sd = NA) # se values

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

# Revel factors here
model_diag_long <- model_diag_long %>% 
  mutate(model = factor(model,
                        levels = c("None", "a", "y", "c",
                                   "a_y", "a_c", "y_a", "y_c", "y_a_c")),
         parameters = factor(parameters,
                             levels = c("rho_a", "rho_y", "rho_c", "log_sigma2"),
                             labels = c( bquote(rho[a]), bquote(rho[y]), bquote(rho[c]), bquote(sigma^2) )))


# Plot parameter estaimtes and WAA RE -------------------------------------------------

# Visualize parameters
(par_plot <- ggplot(model_diag_long, aes(x = model, y = mle_val,
                            ymin = lwr_95, ymax = upr_95,
                            shape = factor(nlminb_conv),
                            color = model)) +
  geom_pointrange(size = 0.95) +
  ggsci::scale_color_jco() +
  scale_shape_manual(values = c(19, 1), 
                     labels = c("Converged", "Not Converged")) +
  facet_wrap(~parameters, scales = "free_y", labeller = label_parsed) +
  guides(color="none") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17)) +
  labs(x = "Model", y = "Parameter Estiamte", shape = "Convergence"))

# Visualize AIC across models
(aic_plot <- ggplot(model_diag_long, aes(x = model, y = AIC, group = 1)) +
  geom_point(aes(color = model),size = 5) +
  geom_line(lty = 2) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17))  +
  labs(x = "Model", y = "AIC"))

png(here("figs", "Par_est_AIC.png"), width = 1300, height = 500)
plot_grid(par_plot, aic_plot, rel_widths = c(0.65, 0.35), 
          align = "hv", axis = "bl")
dev.off()


# Visualize WAA_re --------------------------------------------------------

# Empty dataframe to store values in
WAA_re_df_all <- data.frame()

# extract waa_re values
for(n_fact in 1:nrow(map_factorial)) {
  
  # Extract random parameter values
  rand_par_vals <- models[[n_fact]]$sd_rep$par.random
  
  # coerce these values into a matrix
  WAA_re <- matrix(
      rand_par_vals[names(rand_par_vals) == "WAA_re"],
      ncol = length(ages), nrow = length(years)
    )
  
  # Next, melt into a dataframe
  WAA_re_df <- reshape2::melt(WAA_re)
  colnames(WAA_re_df) <- c("yrs", "ages", "vals")
  WAA_re_df$model <- model_names$model[n_fact] # input model name for indexing
  
  # Now rbind to everything else
  WAA_re_df_all <- rbind(WAA_re_df_all, WAA_re_df)
  
} # end n_fact loop

# Replace blank model with none
WAA_re_df_all$model[WAA_re_df_all$model == ""] <- "None"

# relevel factors
WAA_re_df_all <- WAA_re_df_all %>% 
  mutate(model = factor(model,
                      levels = c("None", "a", "y", "c",
                                 "a_y", "a_c", "y_a", "y_c", "y_a_c")))

# mean value calculation
mean_wt_re <- mean(WAA_re_df_all$vals)

# Now, plot!
plot_list <- list()
pdf(here("figs", "ebs_pollock_WAA_models.pdf"), width = 15, height = 15)

for(i in 1:length(unique(WAA_re_df_all$model))) {
  
  # Iteratively plot these out
    plot_list[[i]] <- ggplot(WAA_re_df_all %>% 
             filter(model == unique(WAA_re_df_all$model)[i]), 
           aes(x = factor(ages), y = factor(yrs), fill = vals)) +
      geom_tile(alpha = 0.9) +
      scale_x_discrete(breaks = seq(3, 15, 3)) +
      scale_y_discrete(breaks = seq(1, 31, 5)) +
      geom_text(aes(label=round(vals,2)), size = 5) +
      scale_fill_gradient2(midpoint = mean_wt_re ) +
      facet_wrap(~model, ncol = 2) +
      theme_bw() +
      labs(x = "Age", y = "Year", fill = "Weight") +
      theme(axis.title = element_text(size = 17),
            axis.text = element_text(size = 15),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            legend.position = "top",
            legend.key.size = unit(0.75, "cm"))
}

plot_list

dev.off()
