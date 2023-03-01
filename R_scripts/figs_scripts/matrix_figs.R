# Creator: Matthew LH. Cheng
# Date 1/8/23
# Purpose: To visualize a triple separable function via random draws using different
# combinations of partial correlations

source(here("R_scripts", "make_precision", "Construct_precision.R"))

library(here)
library(tidyverse)
library(cowplot)

# Plot matrix with different combinastions of correlations -----------------------------------------

# Set up variables here
n_a = 8
n_t = 25

pcorr_age = c(0.1, 0.1, 0.8)
pcorr_year = c(0.1, 0.8, 0.1)
pcorr_cohort = c(0.8, 0.1, 0.1)
var_value = 1e-03

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
  Q = make_precision(n_a, n_t, pcorr_age[i], pcorr_year[i], pcorr_cohort[i], var_value,
                     Var_Type = "Marginal")
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
 

# Random Draws ------------------------------------------------------------

pdf(here("figs", "fig2_rmvnorm_panel.pdf"), height = 6, width = 18)
# Plot!

print(
  mat_plot <- ggplot(corr_all, aes(factor(Year), factor(Age), fill = value )) +
    geom_tile(alpha = 1) +
    facet_wrap(~type, ncol = 3, labeller = label_parsed) +
    scale_fill_gradient2(midpoint = 0, 
                         high = scales::muted("red"),
                         low = scales::muted("blue")) +
    scale_y_discrete(breaks = seq(0, 8, 2)) +
    scale_x_discrete(breaks = seq(0, 25, 5)) +
    labs(x = "Year", y = "Age", fill = "Value") +
    theme_bw() +
    theme(legend.position = "right",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.key.width = unit(0.5, "cm"),
          panel.spacing = unit(3, "lines"), # facet wrap spacing
          strip.text = element_text(size = 17),
          axis.title = element_text(size = 17, vjust = -10),
          axis.text = element_text(size = 15, color = "black"),
          panel.spacing.x = unit(0.5, "lines")) 
)

dev.off()


# Correlation Matrix ------------------------------------------------------

# Turn correlations off for certain partial correlations
pcorr_age = c(0, 0, 0.8)
pcorr_year = c(0, 0.8, 0)
pcorr_cohort = c(0.8, 0, 0)
var_value = 1e-03

# Reducing the dimensions of this to make it more tractable to plot
n_a = 3
n_t = 4

# Get label types
type = c("pa = 0 py = 0 pc = 0.8",
         "pa = 0 py = 0.8 pc = 0",
         "pa = 0.8 pc = 0 pc = 0")

# Make labels bquote to plot greek letters
pc_high = paste(bquote("~rho[a] == 0~~"), bquote("~rho[y] == 0~~"), bquote("~rho[c] == 0.8" )  )
py_high = paste(bquote("~rho[a] == 0~~"), bquote("~rho[y] == 0.8~~"), bquote("~rho[c] == 0" )  )
pa_high = paste(bquote("~rho[a] == 0.8~~"), bquote("~rho[y] == 0~~"), bquote("~rho[c] == 0" )  )

# Loop through to create Y_at matrix dataframe for different corr values
corr_mat_all <- data.frame()
for(i in 1:length(pcorr_age)) {
  
  # Make precision matrix
  Q = make_precision(n_a, n_t, pcorr_age[i], pcorr_year[i], pcorr_cohort[i], var_value,
                     Var_Type = "Marginal")
  # Get dense covariance matrix
  V = solve( Q )
  Vdense = as.matrix(V)
  
  # Simulate random draws here
  corr_mat <- cov2cor(Vdense)
  cor_plot <- corrplot::corrplot(corr_mat, type = 'lower') # correlation plot, but extractin corr vals
  corr_df_vals <- cor_plot$corrPos # get correlatino positions
  corr_df_vals$type <- type[i] # get correlation types
  
  # Get x/year index
  x_idx <- unique(corr_df_vals$x) 
  age_idx <- rep(seq(1, n_a), times = n_t)
  year_idx <- rep(seq(1, n_t), each = n_a)
  x_idx_names <- data.frame(x = x_idx, x_year_idx = year_idx,
                           x_age_idx = age_idx)
  
  # Get y/year index
  y_idx <- unique(corr_df_vals$y)
  age_idx <- rep(seq(1, n_a), times = n_t)
  yr_idx <- rep(seq(1, n_t), each = n_a)
  y_idx_names <- data.frame(y = y_idx, y_age_idx = age_idx, y_year_idx = yr_idx)
  
  # Create age year combinations to differentiate
  corr_df_vals <- corr_df_vals %>% 
    left_join(y_idx_names, by = "y") %>% 
    left_join(x_idx_names, by = "x") %>% 
    mutate(y_idx_names = paste("a", y_age_idx, "y", y_year_idx, sep = ""),
           x_idx_names = paste("a", x_age_idx, "y", x_year_idx, sep = ""))
  
  # now bind into one df
  corr_mat_all <- rbind(corr_mat_all, corr_df_vals)
  
} # end i loop

# Get greek letters 
corr_mat_all <- corr_mat_all %>% 
  mutate(type = factor(type, levels = c("pa = 0.8 pc = 0 pc = 0",
                                        "pa = 0 py = 0.8 pc = 0",
                                        "pa = 0 py = 0 pc = 0.8"),
                       labels = c(pa_high, py_high, pc_high)))

# Relevel factors
corr_mat_all <- corr_mat_all %>% 
  mutate(x_idx_names = factor(x_idx_names,
                              levels = c(paste("a", rep(seq(1, n_a), times = n_t), 
                                               "y", rep(seq(1, n_t), each = n_a), sep = ""))),
         y_idx_names = factor(y_idx_names, 
                              levels = paste("a", rep(seq(n_a, 1), times = n_t), 
                                             "y", rep(seq(n_t, 1), each = n_a), sep = "")))


# Plot out lower triangle of correlation matrix here!
pdf(here("figs", "fig3_corr.pdf"), height = 6, width = 18)
(cor_plot <- ggplot(corr_mat_all, aes(x = factor(x_idx_names), y = factor(y_idx_names), fill = corr,
                         size = corr)) +
  geom_point(alpha = 1, pch = 21) +
  scale_fill_gradient2(low = "white", high = scales::muted("red"),limit = c(0,1)) +
  facet_wrap(~type, labeller = label_parsed, scales = "free_y") +
  labs(x = "", y = "", fill = "Correlation", size = "Correlation") +
  guides(fill = guide_legend(), size = guide_legend())+
  scale_size(range = c(3, 10)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.key.width = unit(1.5, "cm"),
        panel.spacing = unit(3, "lines"), # facet wrap spacing
        strip.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        axis.text.y = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", angle = 90, size = 15),
        panel.spacing.x = unit(0.5, "lines")))
dev.off()        



# Visualize Variance Forms ------------------------------------------------

# Set up variables here
n_a = 3
n_t = 4

pcorr_age = c(0.1, 0.1, 0.8)
pcorr_year = c(0.1, 0.8, 0.1)
pcorr_cohort = c(0.8, 0.1, 0.1)
var_value = 0.01

# Get label types
type = c("pa = 0.1 py = 0.1 pc = 0.8",
         "pa = 0.1 py = 0.8 pc = 0.1",
         "pa = 0.8 pc = 0.1 pc = 0.1")

# Make labels bquote to plot greek letters
pc_high = paste(bquote("~rho[a] == 0.1~~"), bquote("~rho[y] == 0.1~~"), bquote("~rho[c] == 0.8" )  )
py_high = paste(bquote("~rho[a] == 0.1~~"), bquote("~rho[y] == 0.8~~"), bquote("~rho[c] == 0.1" )  )
pa_high = paste(bquote("~rho[a] == 0.8~~"), bquote("~rho[y] == 0.1~~"), bquote("~rho[c] == 0.1" )  )

# Variance Type
var_type <- c("Conditional", "Marginal")

# Store values
Cov_all <- data.frame()

# Loop through to create Y_at matrix dataframe for different corr values
set.seed(666)
for(j in 1:length(var_type)) {
  for(i in 1:length(pcorr_age)) {
    
    # Make precision matrix
    Q = make_precision(n_a, n_t, pcorr_age[i], pcorr_year[i], pcorr_cohort[i], var_value,
                       Var_Type = var_type[j])
    # Get dense covariance matrix
    V = solve( Q )
    Vdense = as.matrix(V)
    
    # Munge dataframe
    Cov_Mat <- reshape::melt(Vdense) %>% 
      rename(Col = X1, Row = X2) %>% 
      mutate(type = type[i],
             var_type = var_type[j])
    
    # now bind into one df
    Cov_all <- rbind(Cov_Mat, Cov_all)
    
  } # end i loop
} # end j loop

# Get greek letters 
Cov_all <- Cov_all %>% 
  mutate(type = factor(type, levels = c("pa = 0.8 pc = 0.1 pc = 0.1",
                                        "pa = 0.1 py = 0.8 pc = 0.1",
                                        "pa = 0.1 py = 0.1 pc = 0.8"),
                       labels = c(pa_high, py_high, pc_high)),
         diag = ifelse(Col == Row, "Variance", "Covariance")) # diagonal matrix


pdf(here("figs", "figc1_varcov.pdf"), height = 15, width = 25)
# Now plot this out!
ggplot(Cov_all, aes(x = factor(Col), y = factor(Row), fill = value)) +
  geom_tile(alpha = 1) +
  scale_fill_gradient2(high = scales::muted("red")) +
  geom_text(aes(label = round(value, 2), color = diag), size = 6.5) +
  facet_grid(var_type ~ type, labeller = label_parsed) +
  labs(x = "Columns", y = "Rows", fill = "Values", color = "Matrix Component") +
  scale_color_manual(values = c("black", "blue")) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 21),
        legend.key.width = unit(0.5, "cm"),
        panel.spacing = unit(0.1, "lines"), # facet wrap spacing
        strip.text = element_text(size = 21),
        axis.title = element_text(size = 21),
        axis.text.y = element_text(color = "black", size = 19),
        axis.text.x = element_text(color = "black", angle = 90, size = 19))
dev.off() 
