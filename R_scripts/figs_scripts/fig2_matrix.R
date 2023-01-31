# Creator: Matthew LH. Cheng
# Date 1/8/23
# Purpose: To visualize a triple separable function via random draws using different
# combinations of partial correlations

source(here("R_scripts", "make_precision", "Construct_precision_2023-01-02.R"))

library(here)
library(tidyverse)

# Plot matrix with different combinastions of correlations -----------------------------------------

# Set up variables here
n_a = 8
n_t = 25

pcorr_age = c(0.1, 0.1, 0.8)
pcorr_year = c(0.1, 0.7, 0.1)
pcorr_cohort = c(0.8, 0.1, 0.1)
margvar = 1e-03

# Get label types
type = c("pa = 0.1 py = 0.1 pc = 0.8",
         "pa = 0.1 py = 0.8 pc = 0.1",
         "pa = 0.8 pc = 0.1 pc = 0.1")

# Make labels bquote to plot greek letters
pc_high = paste(bquote("~rho[a] == 0.1~~"), bquote("~rho[y] == 0.1~~"), bquote("~rho[c] == 0.8" )  )
py_high = paste(bquote("~rho[a] == 0.1~~"), bquote("~rho[y] == 0.8~~"), bquote("~rho[c] == 0.1" )  )
pa_high = paste(bquote("~rho[a] == 0.8~~"), bquote("~rho[y] == 0.1~~"), bquote("~rho[c] == 0.1" )  )

# Store values
corr_all <- data.frame()

# Loop through to create Y_at matrix dataframe for different corr values
set.seed(666)
for(i in 1:length(pcorr_age)) {
  
  # Make precision matrix
  Q = make_precision(n_a, n_t, pcorr_age[i], pcorr_year[i], pcorr_cohort[i], margvar)
  # Get dense covariance matrix
  V = solve( Q )
  Vdense = as.matrix(V)
  
  # Simulate random draws here
  Y_at = matrix( rmvnorm(n=1, mean=rep(0,n_a*n_t), sigma=Vdense), nrow=n_a, ncol=n_t )
  
  # Munge dataframe
  Y_mat <- reshape::melt(Y_at) %>% 
    rename(Age = X1, Year = X2) %>% 
    mutate(type = type[i])
  
  # now bind into one df
  corr_all <- rbind(Y_mat, corr_all)
  
} # end i loop

# Get greek letters 
corr_all <- corr_all %>% 
  mutate(type = factor(type, levels = c("pa = 0.8 pc = 0.1 pc = 0.1",
                                        "pa = 0.1 py = 0.8 pc = 0.1",
                                        "pa = 0.1 py = 0.1 pc = 0.8"),
                       labels = c(pa_high, py_high, pc_high)))

pdf(here("figs", "fig2_rmvnorm_panel.pdf"), height = 5, width = 18)
# Plot!
print(
  ggplot(corr_all, aes(factor(Year), factor(Age), fill = value )) +
    geom_tile(alpha = 1.5) +
    facet_wrap(~type, ncol = 3, labeller = label_parsed) +
    scale_fill_gradient2() +
    scale_y_discrete(breaks = seq(0, 8, 2)) +
    scale_x_discrete(breaks = seq(0, 25, 5)) +
    labs(x = "Year", y = "Age", fill = "Value") +
    theme_bw() +
    theme(legend.position = "top",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.key.width = unit(1.5, "cm"),
          panel.spacing = unit(3, "lines"), # facet wrap spacing
          strip.text = element_text(size = 17),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black")) 
)
dev.off()
