# Purpose: To implement and test triple separability .cpp file
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date: 1.15.23

# set up ------------------------------------------------------------------
setwd(R'(C:\Users\James.Thorson\Desktop\Git\Triple_Separability)')

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
Xse_at <- t(as.matrix(waa_std_df[,-1])) # removing first col (year column)

# Create an index for ages and years to feed into TMB, which helps construct the precision matrix
ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), 
                                  "year" = seq_len(length(years)) ))


# Set up TMB Model --------------------------------------------------------

# Now, input these components into a data list
data <- list( years = years,
              ages = ages,
              X_at = X_at,
              Xse_at = Xse_at,
              ay_Index = ay_Index )

# Input parameters into a list
parameters <- list( rho_y = rbeta(1, 1, 1),
                    rho_a = rbeta(1, 1, 1),
                    rho_c = rbeta(1, 1, 1),
                    log_sigma2 = rbeta(1, 1, 1),
                    ln_L0 = log(0.1),
                    ln_Linf = log(1),  # Change to reasonable value, only a*Linf is identified using WAA data
                    ln_k = log(0.2),
                    ln_alpha = log(1),
                    ln_beta = log(3),   # Fix at isometric
                    Y_at = array(0,dim=dim(X_at)) )

#
map = list( "ln_Linf" = factor(NA),
            "ln_beta" = factor(NA) )

# Now, make AD model function
waa_model <- MakeADFun(data = data, parameters = parameters, 
                       random = "Y_at",
                       DLL = "triple_sep_waa",
                       map = map)
report = waa_model$report()
report$mu_at

# Now, optimize the function
waa_optim <- stats::nlminb(waa_model$par, waa_model$fn, waa_model$gr,  
                           control = list(iter.max = 1e5, eval.max = 1e5))
report = waa_model$report()
report$mu_at

# Get report
waa_model$rep <- waa_model$report(waa_model$env$last.par.best)

# Get sd report
sd_rep <- sdreport(waa_model)

# Visualize sparse matrix
Matrix::image(waa_model$env$spHess(random=TRUE))

# Check convergence
waa_optim$convergence == 0 
sd_rep$pdHess == TRUE  
max(abs(sd_rep$gradient.fixed))
waa_model$rep$jnLL # joint nLL seems very large...
waa_optim$objective 

# Plot covariance
Q = waa_model$report()$Q
V = solve(Q)
R = cov2cor(V)
P_at = matrix( R[,14], nrow=length(ages), ncol=length(years) )
image(t(P_at))

# Extract values ----------------------------------------------------------

# Extract WAA random effects
#WAA_re <- matrix(
#  sd_rep$par.random[names(sd_rep$par.random) == "WAA_re"],
#  ncol = length(ages), nrow = length(years)
#)
WAA_re <- t(waa_model$env$parList()$Y_at)

# Munge WAA for plotting
WAA_re <- reshape2::melt(WAA_re) 
colnames(WAA_re) <- c("yrs", "ages", "vals")
WAA_re <- WAA_re %>% 
  mutate(Type = "Random")

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


# Visualize! --------------------------------------------------------------

ggplot(WAA_all, aes(x = factor(ages), y = factor(yrs), fill = vals)) +
  geom_tile(alpha = 0.9) +
  scale_x_discrete(breaks = seq(3, 15, 3)) +
  scale_y_discrete(breaks = seq(1, 31, 5)) +
  geom_text(aes(label=round(vals,2)), size = 5) +
  scale_fill_gradient2(midpoint = mean_wt_re ) +
  facet_wrap(~Type) +
  theme_bw() +
  labs(x = "Age", y = "Year") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))
