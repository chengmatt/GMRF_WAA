# Purpose: To visualize GMRF models for walleye Pollock (for conditional variance)
# Creator: Matthew LH. Cheng
# Date 1/19/23


# Set up -----------------------------------------------------------------

# load in packages
library(here)
library(tidyverse)
library(cowplot)
library(Matrix)

# Load in models
load(here("output", "cond_var_ebs_pollock_waa_models.RData")) # obj is called models

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
   geom_pointrange(size = 1.1) +
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
   labs(x = "Model", y = "Parameter Estiamte", shape = "Convergence"))

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
plot_grid(par_plot, aic_plot, rel_widths = c(0.5, 0.5), 
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
  rand_par_vals <- exp(models[[n_fact]]$sd_rep$par.random)
  
  # coerce these values into a matrix
  WAA_re <- matrix(
    t(exp(models[[n_fact]]$env$parList()$ln_Y_at)),
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
                                   "a_y", "a_c", "y_a", "y_c", "y_a_c"),
                        labels = c("None", "Age", "Year", "Cohort", 
                                   "Age+Year", "Age+Cohort", "Year+Age", "Year+Cohort", "Year+Age+Cohort")))

# Relevel for mean WAA df
mean_WAA_all <- mean_WAA_all %>% 
  mutate(model = factor(model,
                        levels = c("None", "a", "y", "c",
                                   "a_y", "a_c", "y_a", "y_c", "y_a_c"),
                        labels = c("None", "Age", "Year", "Cohort", 
                                   "Age+Year", "Age+Cohort", "Year+Age", "Year+Cohort", "Year+Age+Cohort")))

# Compute anomaly relative to the mean
WAA_re_df_all <- WAA_re_df_all %>% 
  left_join(mean_WAA_all, by = c("ages", "model")) %>% 
  mutate(anom = (vals - mean_waa) / mean_waa)

tile_plot <- ggplot(WAA_re_df_all, 
                    aes(y = factor(ages), x = factor(yrs), fill = anom)) +
  geom_tile(alpha = 1) +
  scale_y_discrete(breaks = seq(3, 15, 3)) +
  scale_x_discrete(breaks = seq(1990, 2020, 5)) +
  scale_fill_gradient2(midpoint = mean(WAA_re_df_all$anom), 
                       breaks = seq(-0.1, 0.6, 0.2),
                       high = scales::muted("red"), 
                       low = scales::muted("blue")) +
  theme_bw() +
  facet_wrap(~model, ncol = 4) +
  labs(x = "Year", y = "Age", fill = "Anomaly relative to mean WAA") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 17),
        legend.position = "top",
        legend.key.width = unit(1.5, "cm"))

line_plot <- ggplot(WAA_re_df_all %>% 
                      filter(ages %in% c(seq(3, 15, 2))), 
                      aes(x = factor(yrs), y = vals, color = factor(ages),
                      group = factor(ages))) +
  geom_line(alpha = 1, size = 1.6) +
  # geom_text(aes(label = ages), size = 4.5) +
  scale_x_discrete(breaks = seq(1990, 2020, 5)) +
  ggsci::scale_color_jco( ) +
  theme_bw() +
  facet_wrap(~model, ncol = 4) +
  guides(color=guide_legend(ncol=3)) +
  labs(x = "Year", y = "Weight", color = "Ages") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 17),
        legend.position = c(0.055, 0.92),
        legend.background = element_blank(),
        legend.key.width = unit(0.75, "cm"))

# Now, plot!
png(here("figs", "ebs_pollock_WAA_models_tile.png"), width = 1400, height = 1000)

plot_grid(tile_plot, line_plot, rel_heights = c(0.8, 1), axis = "bl", align = "hv",
          ncol = 1, labels = c("A", "B"), hjust = c(-2, -2),  
          vjust = c(3, 0), label_size = 23)

dev.off()


# Uncertainty -------------------------------------------------------------

# Get CV
mod_none_cv <- (models[[1]]$sd_rep$sd / exp(models[[1]]$sd_rep$value)) * 100
mod_all_cv <- (models[[8]]$sd_rep$sd / exp(models[[8]]$sd_rep$value)) * 100

# Coerce into dataframe
mod_none_df <- reshape2::melt(matrix(mod_none_cv, 
                              nrow = length(years), ncol = length(ages))) %>% 
  mutate(Model = "None")

mod_all_df <- reshape2::melt(matrix(mod_all_cv, 
                         nrow = length(years), ncol = length(ages))) %>% 
  mutate(Model = "Year+Age+Cohort")

# Bind to df
mod_cv_df <- rbind(mod_none_df, mod_all_df)
colnames(mod_cv_df) <- c("Year", "Age", "Vals", "Model")

# Now plot the differences!
mod_cv_df_wide <- mod_cv_df %>% 
  pivot_wider(names_from = "Model", values_from = "Vals") %>% 
  mutate(Difference = None - `Year+Age+Cohort`)
  
pdf(here("figs", "cv_waa_tile.pdf"), width = 15, height = 10)
ggplot(mod_cv_df_wide, aes(x = factor(Year + 1990) , y = factor(Age + 2), 
                           fill = Difference )) +
  geom_tile() +
  geom_text(aes(label = round(Difference, 2))) +
  scale_y_discrete(breaks = seq(3, 15, 3)) +
  scale_x_discrete(breaks = seq(1990, 2020, 5)) +
  scale_fill_gradient2(midpoint = 0, 
                       high = scales::muted("red"), 
                       low = scales::muted("blue")) +
  theme_bw() +
  labs(x = "Year", y = "Age", 
       fill = "CV(%) difference (Model None - Year+Age+Cohort)") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 17),
        legend.position = "top",
        legend.key.width = unit(1.5, "cm"))
dev.off()
