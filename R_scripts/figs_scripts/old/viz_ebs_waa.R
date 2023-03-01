# Purpose: Visualize the EBS pollock WAA matrix
# Creator: Matthew LH. Cheng
# Date 1/15/22


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)

# Load in data
ebs_df <- read.csv(here("data", "ebs_waa.csv"))
ebs_std <- read.csv(here("data", "ebs_waa_std.csv"))

# Plot! -------------------------------------------------------------------

# Get mean weight
mean_wt <- ebs_df %>% 
  pivot_longer(!c(year, source), values_to = "wt", names_to = "ages") %>% 
  mutate(ages = parse_number(ages)) %>% 
  filter(source == "fishery") %>% 
  summarize(mean = mean(wt))

# get mean std
std_wt <- ebs_std %>% 
  pivot_longer(!c(year, source), values_to = "std", names_to = "ages") %>% 
  mutate(ages = parse_number(ages)) %>% 
  filter(source == "fishery") %>% 
  summarize(mean = mean(std))

png(filename = here("figs", "ebs_pollock_waa_mat.png"), width = 800, height = 1000)

ebs_df %>% 
  pivot_longer(!c(year, source), values_to = "wt", names_to = "ages") %>% 
  mutate(ages = parse_number(ages)) %>% 
  filter(source == "fishery") %>% 
  ggplot(aes(x = factor(ages), y = factor(year), fill = wt)) +
  geom_tile(alpha = 0.9) +
  scale_x_discrete(breaks = seq(3, 15, 3)) +
  scale_y_discrete(breaks = seq(1991, 2021, 5)) +
  geom_text(aes(label=round(wt,2)), size = 5) + 
  scale_fill_gradient2(midpoint = mean_wt$mean) +
  theme_bw() +
  labs(x = "Age", y = "Year") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))

dev.off()


# Visualize stds
png(filename = here("figs", "ebs_pollock_waa__std_mat.png"), width = 800, height = 1000)
ebs_std %>% 
  pivot_longer(!c(year, source), values_to = "std", names_to = "ages") %>% 
  mutate(ages = parse_number(ages)) %>% 
  filter(source == "fishery") %>% 
  ggplot(aes(x = factor(ages), y = factor(year), fill = std)) +
  geom_tile(alpha = 0.9) +
  scale_x_discrete(breaks = seq(3, 15, 3)) +
  scale_y_discrete(breaks = seq(1991, 2021, 5)) +
  geom_text(aes(label=round(std,2)), size = 5) + 
  scale_fill_gradient2(midpoint = 0.3) +
  theme_bw() +
  labs(x = "Age", y = "Year") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))
dev.off()

# Visualize as line plot
ebs_df %>% 
  pivot_longer(!c(year, source), values_to = "wt", names_to = "ages") %>% 
  mutate(ages = parse_number(ages) + 2) %>% 
  filter(source == "fishery",
         ages %in% c(seq(3, 15, 2))) %>% 
  ggplot(aes(x = year, y = wt, color = factor(ages))) +
  geom_line(alpha = 1, size = 2) +
  ggsci::scale_color_jco( ) +
  guides(color=guide_legend(ncol=3, override.aes = list(size = 5))) +
  scale_fill_gradient2(midpoint = mean_wt$mean) +
  theme_bw() +
  labs(x = "Year", y = "Weight", color = "Ages") +
  theme(axis.title = element_text(size = 23),
        axis.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 21),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 17),
        legend.position = c(0.08, 0.92),
        legend.background = element_blank(),
        legend.key.width = unit(0.75, "cm"))
