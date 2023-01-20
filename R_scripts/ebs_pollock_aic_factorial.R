# Purpose: To run a series of triple separable models testing for evidence of
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
             ay_Index = ay_Index, Var_Param = 0) # Conditional Variance

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

save(models, file = here("output", "ebs_pollock_waa_models.RData"))


# Visualize WAA_re --------------------------------------------------------

# Empty dataframe to store values in
WAA_re_df_all <- data.frame()

# extract waa_re values
for(n_fact in 1:nrow(map_factorial)) {
  
  # Extract random parameter values
  rand_par_vals <- models[[n_fact]]$sd_rep$par.random
  
  # coerce these values into a matrix
  WAA_re <- matrix(
      rand_par_vals[names(rand_par_vals) == "Y_at"],
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

# Calculate anomaly relative to the mean weight-at-age
mean_waa <- reshape2::melt(models[[1]]$rep$mu_at[,1]) %>% 
  mutate(ages = 1:13) %>% 
  rename(mean_waa = value)

# Compute anomaly relative to the mean
WAA_re_df_all <- WAA_re_df_all %>% 
  left_join(mean_waa, by = "ages") %>% 
  mutate(anom = (vals - mean_waa) / mean_waa)

# Now, plot!
pdf(here("figs", "ebs_pollock_WAA_models_tile.pdf"), width = 15, height = 15)

 tile_plot <- ggplot(WAA_re_df_all, 
         aes(x = factor(ages), y = factor(yrs), fill = anom)) +
    geom_tile(alpha = 0.9) +
    scale_x_discrete(breaks = seq(3, 15, 3)) +
    scale_y_discrete(breaks = seq(1, 31, 5)) +
    scale_fill_gradient2( ) +
    theme_bw() +
    facet_wrap(~model, ncol = 4) +
    labs(x = "Age", y = "Year", fill = "Anomaly relative to mean WAA") +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 17),
          legend.position = "top",
          legend.key.width = unit(1.5, "cm"))
 
 line_plot <- ggplot(WAA_re_df_all %>% 
          filter(ages %in% c(seq(1, 13, 3))), 
        aes(x = factor(yrs), y = vals, color = factor(ages),
            group = factor(ages))) +
   geom_line(alpha = 0.75, size = 1.3) +
   geom_text(aes(label = ages), size = 5) +
   scale_x_discrete(breaks = seq(1, 31, 5)) +
   ggsci::scale_color_jco( ) +
   theme_bw() +
   facet_wrap(~model, ncol = 4) +
   labs(x = "Year", y = "Wieght") +
   theme(axis.title = element_text(size = 17),
         axis.text = element_text(size = 15),
         legend.title = element_text(size = 17),
         legend.text = element_text(size = 15),
         strip.text = element_text(size = 17),
         legend.position = "none",
         legend.key.width = unit(1.5, "cm"))
 
 plot_grid(tile_plot, line_plot, rel_heights = c(1, 1),
           ncol = 1)


dev.off()
