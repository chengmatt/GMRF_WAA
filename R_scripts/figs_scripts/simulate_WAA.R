# Creator: Matthew LH. Cheng
# Purpose: To simulate WAA with age, cohort, and year effects
# Date 1/8/23

# packages
library(here)
library(tidyverse)
library(mvtnorm) 

# Load in precision matrix constructor function
source(here("scripts", "Construct_precision_2023-01-02.R"))

# Set up variables here
n_a = 15 # ages
n_t = 35 # years

pcorr_age = 0.1 # partial correlations along age
pcorr_year = 0.1 # partial correlations along year
pcorr_cohort = 0.8 # partial correlations along cohort
margvar = 0.3 # marginal variance (diagonals)


# Set up precision matrix -------------------------------------------------

# We need a vector of mean weights with length (n = age x year)
ages = 1:n_a

# Arbrtirary parameter values for now
winf = 10 
k = 1.5
allom_fct = 3

# Calculate mean weight
mean_weight = winf * (1 - exp(-k * ages)) ^ allom_fct # mean ages 

# Convert the above into a mean vector of ages dimensioned by age x year x 1)
weight_vec = as.vector(rep(mean_weight, times = n_t))

# Make precision matrix
sparse_Q = make_precision(n_a = n_a, n_t = n_t, pcorr_age = pcorr_age,
                          pcorr_year = pcorr_year, pcorr_cohort = pcorr_cohort,
                          margvar = margvar, what = "Q")

# Now solve to make our covariance matrix
dense_cov_mat <- as.matrix(solve(sparse_Q))

# Now, simulate random draws with this vector
weight_a_t = matrix(rmvnorm(n = 1, mean = weight_vec, sigma = dense_cov_mat),
                    nrow=n_a, ncol=n_t)

# Each of these columns represent a single year of mean stock weights


# Plot simulated data -----------------------------------------------------

wt_df = reshape2::melt(weight_a_t) %>% 
  rename(Ages = Var1, Years = Var2, Weights = value)

png(here("figs", "fig3_placeholder_waa_lines.png"), height = 600, width = 750)

# Weight-at-age plot across years
ggplot(wt_df, aes(x = Ages, y = Weights, group = Years, color = Years)) +
  geom_line(alpha = 0.65) +
  # Mean weight vector
  geom_line(aes(x = Ages, y = weight_vec), color = "black", size = 1.25, lty = 2) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(x = "Ages", y = "Weights")

dev.off()

png(here("figs", "fig3_placeholder_waa_heat.png"), height = 600, width = 750)

# Weight-at-age tile plot
ggplot(wt_df, aes(factor(Years), factor(Ages), fill = Weights)) +
  geom_tile() +
  scale_y_discrete(breaks = seq(0, 15, 3)) +
  scale_x_discrete(breaks = seq(0, 35, 5)) +
  scale_fill_gradient2(midpoint = 7.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(x = "Years", y = "Ages")

dev.off()
