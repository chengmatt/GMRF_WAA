# Purpose: To visualize GMRF models for walleye Pollock
# Creator: Matthew LH. Cheng
# Date 1/19/23


# Set up -----------------------------------------------------------------

# load in packages
library(here)
library(tidyverse)
library(cowplot)

# Load in models
load(here("output", "ebs_pollock_waa_models.RData")) # obj is called models

# Load in model diagnostics csv
model_diag <- read.csv(here("output", "model_diag_vals.csv"))

# Load in WAA matrix (only use fishery data)
waa_df <- read.csv(here("data", "ebs_waa.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Load in std for WAA matrix
waa_std_df <- read.csv(here("data", "ebs_waa_std.csv")) %>% 
  filter(source == "fishery") %>% 
  dplyr::select(-source)

# Generate full factorial design for naming models
map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>% 
  data.frame()

# Relevel factors here
model_diag <- model_diag %>% 
  mutate(model = factor(model, levels = c("None", "a", "y", "c",
                                          "a_y", "a_c", "y_a", "y_c", "y_a_c")),
         parameters = factor(parameters,
                             levels = c("rho_a", "rho_y", 
                                        "rho_c", "log_sigma2"),
                             labels = c( bquote(rho[a]), bquote(rho[y]), 
                                         bquote(rho[c]), bquote(sigma^2) )))

# Create model names to differentiate models
model_names <- map_factorial %>% 
  mutate(rho_y_lab = case_when(rho_y == 1 ~ "y"),
         rho_a_lab = case_when(rho_a == 1 ~ "a"),
         rho_c_lab = case_when(rho_c == 1 ~ "c")) %>% 
  dplyr::select(rho_y_lab, rho_a_lab, rho_c_lab) %>% 
  dplyr::rowwise() %>% 
  tidyr::unite('model', na.rm = TRUE)

# Dimensions
years <- 1991:2021
ages <- 3:15

# Plot parameter estimates and WAA RE -------------------------------------------------

# Visualize parameters
(par_plot <- ggplot(model_diag, aes(x = model, y = mle_val,
                                         ymin = lwr_95, ymax = upr_95,
                                         shape = factor(nlminb_conv),
                                         color = model)) +
   geom_pointrange(size = 0.95) +
   ggsci::scale_color_jco() +
   scale_shape_manual(values = c(19, 1), 
                      labels = c("Converged", "Not Converged")) +
   facet_wrap(~parameters, scales = "free_y", labeller = label_parsed, ncol = 1) +
   guides(color="none") +
   theme_bw() +
   theme(legend.position = 'none',
         axis.title = element_text(size = 17),
         axis.text = element_text(size = 15, color = "black"),
         strip.text = element_text(size = 17)) +
   labs(x = "", y = "Parameter Estiamte", shape = "Convergence"))

# Visualize AIC across models
(aic_plot <- ggplot(model_diag, aes(x = model, y = AIC, group = 1)) +
    geom_line(lty = 2, size = 1.3) +
    geom_point(aes(color = model),size = 5) +
    ggsci::scale_color_jco() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          strip.text = element_text(size = 17))  +
    labs(x = "Model", y = "AIC"))

pdf(here("figs", "Par_est_AIC.pdf"), width = 17, height = 10)
plot_grid(aic_plot, par_plot, rel_widths = c(0.5, 0.5), 
          align = "hv", axis = "bl", ncol = 2,
          labels = c("A", "B"), label_size = 23, hjust = -0.5)
dev.off()
 

# Visualize WAA_re --------------------------------------------------------

# Empty dataframe to store values in
WAA_re_df_all <- data.frame()
mean_WAA_all <- data.frame()

# extract waa_re values
for(n_fact in 1:nrow(map_factorial)) {
  
  # Extract random parameter values
  rand_par_vals <- models[[n_fact]]$sd_rep$par.random
  
  # coerce these values into a matrix
  WAA_re <- matrix(
    t(models[[n_fact]]$env$parList()$Y_at),
    ncol = length(ages), nrow = length(years)
  )
  
  # Next, melt into a dataframe
  WAA_re_df <- reshape2::melt(WAA_re)
  colnames(WAA_re_df) <- c("yrs", "ages", "vals")
  WAA_re_df$model <- model_names$model[n_fact] # input model name for indexing
  WAA_re_df$ages <- WAA_re_df$ages + 2 # Adding the true start age back
  WAA_re_df$yrs <- WAA_re_df$yrs + 1990 # Adding the true start year back
  
  # Do the same but for mean WAA
  mean_waa <- reshape2::melt(models[[n_fact]]$rep$mu_at[,1]) %>% 
  mutate(ages = 3:15,
         model = model_names$model[n_fact]) %>% 
  rename(mean_waa = value)
  
  
  # Now rbind to everything else
  WAA_re_df_all <- rbind(WAA_re_df_all, WAA_re_df)
  mean_WAA_all <- rbind(mean_waa, mean_WAA_all)
  
} # end n_fact loop

# Replace blank model with none
WAA_re_df_all$model[WAA_re_df_all$model == ""] <- "None"
mean_WAA_all$model[mean_WAA_all$model == ""] <- "None"

# relevel factors
WAA_re_df_all <- WAA_re_df_all %>% 
  mutate(model = factor(model,
                        levels = c("None", "a", "y", "c",
                                   "a_y", "a_c", "y_a", "y_c", "y_a_c")))

# Relevel for mean WAA df
mean_WAA_all <- mean_WAA_all %>% 
  mutate(model = factor(model,
                        levels = c("None", "a", "y", "c",
                                   "a_y", "a_c", "y_a", "y_c", "y_a_c")))

# Compute anomaly relative to the mean
WAA_re_df_all <- WAA_re_df_all %>% 
  left_join(mean_WAA_all, by = c("ages", "model")) %>% 
  mutate(anom = (vals - mean_waa) / mean_waa)

tile_plot <- ggplot(WAA_re_df_all, 
                    aes(x = factor(ages), y = factor(yrs), fill = anom)) +
  geom_tile(alpha = 0.9) +
  scale_x_discrete(breaks = seq(3, 15, 3)) +
  scale_y_discrete(breaks = seq(1990, 2020, 5)) +
  scale_fill_gradient2( ) +
  theme_bw() +
  facet_wrap(~model, ncol = 4) +
  labs(x = "Age", y = "Year", fill = "Anomaly relative to mean WAA") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17),
        legend.position = "top",
        legend.key.width = unit(1.5, "cm"))

line_plot <- ggplot(WAA_re_df_all %>% 
                      filter(ages %in% c(seq(3, 15, 3))), 
                    aes(x = factor(yrs), y = vals, color = factor(ages),
                        group = factor(ages))) +
  geom_line(alpha = 0.75, size = 1.3) +
  geom_text(aes(label = ages), size = 5.5) +
  scale_x_discrete(breaks = seq(1990, 2020, 5)) +
  ggsci::scale_color_jco( ) +
  theme_bw() +
  facet_wrap(~model, ncol = 4) +
  labs(x = "Year", y = "Wieght") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17),
        legend.position = "none",
        legend.key.width = unit(1.5, "cm"))

# Now, plot!
pdf(here("figs", "ebs_pollock_WAA_models_tile.pdf"), width = 18, height = 15)

plot_grid(tile_plot, line_plot, rel_heights = c(1, 1), axis = "bl", align = "hv",
          ncol = 1, labels = c("A", "B"), hjust = c(-2, -2),  
          vjust = c(3, 0), label_size = 23)

dev.off()
   