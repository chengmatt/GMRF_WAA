# Purpose: To implement and test triple separability .cpp file
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date: 1.15.23

# set up ------------------------------------------------------------------
# setwd(R'(C:\Users\James.Thorson\Desktop\Git\Triple_Separability)')

library(here)
library(tidyverse)
library(TMB)

# Load in WAA matrix (only use fishery data)
waa_df <- read.csv(here("data", "ebs_waa.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Load in std for WAA matrix
waa_std_df <- read.csv(here("data", "ebs_waa_std.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Compile and load in model
setwd(here("src"))
compile("triple_sep_waa.cpp")
dyn.load(dynlib("triple_sep_waa"))


# Set up TMB data ----------------------------------------

# Years
years <- waa_df$year

# Ages (goes from age 3 - 15+)
ages <- parse_number(colnames(waa_df)[-1])

# Read in data weight at age matrix
X_at <- t(as.matrix(waa_df[,-1])) # removing first col (year column)

# Read in standard deviations for weight at age matrix
Xse_at <- t(as.matrix(waa_std_df[,c(-1)])) # removing first col (year column) 

# Convert to CV
Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )

# Now convert back to sd in lognormal space
Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))

# Create an index for ages and years to feed into TMB, which helps construct the precision matrix
ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), 
                                  "year" = seq_len(length(years)) ))


# Set up TMB Model --------------------------------------------------------

# Now, input these components into a data list
data <- list( years = years,
              ages = ages,
              X_at = X_at,
              Xsd_at = Xsd_at,
              ay_Index = ay_Index,
              Var_Param = 0) # Var_Param == 0 Conditional, == 1 Marginal

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
                    ln_Y_at = array(0,dim=dim(X_at))) 

# Turn params off
map = list( "ln_Linf" = factor(NA),
            "ln_beta" = factor(NA)
)

compile("triple_sep_waa.cpp")
dyn.load(dynlib("triple_sep_waa"))

# Now, make AD model function
waa_model <- MakeADFun(data = data, parameters = parameters, 
                       random = "ln_Y_at",
                       DLL = "triple_sep_waa",
                       map = map, silent = FALSE)

report = waa_model$report()
diag(solve(report$Q_sparse))

# Now, optimize the function
waa_optim <- stats::nlminb(waa_model$par, waa_model$fn, waa_model$gr,  
                           control = list(iter.max = 1e5, eval.max = 1e5))

report = waa_model$report()
plot(report$mu_at[,1], type = "l")

# Get sd report
sd_rep <- sdreport(waa_model)

# Check marginal variance
diag(solve(report$Q_sparse))

# Visualize sparse matrix
Matrix::image(waa_model$env$spHess(random=TRUE))

# Check convergence
waa_optim$convergence == 0 
sd_rep$pdHess == TRUE  
max(abs(sd_rep$gradient.fixed))
waa_model$report()$jnLL 
waa_optim$objective 


# Plot covariance for model checking --------------------------------------

# Model covariance
Q = waa_model$report()$Q
V = solve(Q)
diag(V)
R = cov2cor(V)
P_at = matrix( R[,403], nrow=length(ages), ncol=length(years) )
image(t(P_at)) 


# Extract values ----------------------------------------------------------

# Extract WAA random effects
WAA_re <- t(waa_model$env$parList()$Y_at)

# Munge WAA for plotting
WAA_re <- reshape2::melt(WAA_re) 
colnames(WAA_re) <- c("yrs", "ages", "vals")
WAA_re <- WAA_re %>% mutate(Type = "Random")

# Get mean predicted weights
mean_wt_re <- mean(WAA_re$vals)

# Get empirical WAA here as well
WAA_emp <- as.matrix(waa_df[,-1])
WAA_emp <- reshape2::melt(WAA_emp) 
colnames(WAA_emp) <- c("yrs", "ages", "vals")

# Now, do some quick munging
WAA_emp <- WAA_emp %>% 
  mutate(ages = parse_number(paste(ages)),
         Type = "Empirical")

WAA_all <- rbind(WAA_emp, WAA_re)


# Visualize for model checking --------------------------------------------------------------

# Calculate anomaly relative to the mean weight-at-age
mean_waa <- reshape2::melt(report$mu_at[,1]) %>% 
  mutate(ages = 1:13) %>% 
  rename(mean_waa = value)

# Compute anomaly relative to the mean
WAA_all <- WAA_all %>% 
  left_join(mean_waa, by = "ages") %>% 
  mutate(anom = (vals - mean_waa) / mean_waa)

ggplot(WAA_all %>% filter(Type == "Random"), 
       aes(x = factor(ages), y = factor(yrs), fill = anom)) +
  geom_tile(alpha = 0.9) +
  scale_x_discrete(breaks = seq(3, 15, 3)) +
  scale_y_discrete(breaks = seq(1, 31, 5)) +
  geom_text(aes(label=round(vals,2)), size = 5) +
  scale_fill_gradient2(midpoint = mean(WAA_all$anom) ) +
  facet_wrap(~Type) +
  theme_bw() +
  labs(x = "Age", y = "Year") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))

ggplot(WAA_all %>% filter(Type == "Random"), 
       aes(x = factor(yrs), y = vals, color = factor(ages),
           group = factor(ages))) +
  geom_text(aes(label=round(ages,2)), size = 4.5) +
  geom_line(alpha = 0.85) +
  theme_bw() +
  labs(x = "Year", y = "Weight") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "none")
