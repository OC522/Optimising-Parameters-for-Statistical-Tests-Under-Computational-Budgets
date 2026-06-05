# Chapter 4 plotting notes.

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(scales)

# Load data.
df <- read.csv("results_main_validonly.csv", stringsAsFactors = FALSE)

# Keep the baseline methods for Chapter 4.
baseline_methods <- c("Uniform", "Random", "SA", "SH")
df <- df %>%
  filter(method %in% baseline_methods)

# Keep plotting order stable.
df <- df %>%
  mutate(
    method = factor(method, levels = baseline_methods),
    theta_set_size = factor(theta_set_size, levels = c(64, 128, 256, 512)),
    B_tune_sec = factor(B_tune_sec, levels = c(40, 80, 160, 320))
  )

# Make the figures folder if needed.
if (!dir.exists("figures")) dir.create("figures")

# Summaries by regime.
summary_regime <- df %>%
  group_by(theta_set_size, B_tune_sec, method) %>%
  summarise(
    mean_power = mean(final_power_hat, na.rm = TRUE),
    se_power = sd(final_power_hat, na.rm = TRUE) / sqrt(sum(!is.na(final_power_hat))),
    
    mean_valid_rate = mean(valid_flag, na.rm = TRUE),
    
    mean_valid_power = mean(valid_power, na.rm = TRUE),
    se_valid_power = sd(valid_power, na.rm = TRUE) / sqrt(sum(!is.na(valid_power))),
    
    mean_type1 = mean(final_type1_hat, na.rm = TRUE),
    mean_c0 = mean(final_c0_mean, na.rm = TRUE),
    
    mean_n_eval = mean(n_theta_eval, na.rm = TRUE),
    .groups = "drop"
  )

# Rank within each regime.
rank_df <- summary_regime %>%
  group_by(theta_set_size, B_tune_sec) %>%
  arrange(desc(mean_valid_power), .by_group = TRUE) %>%
  mutate(rank_valid_power = row_number()) %>%
  ungroup()

# Shared plot theme.
theme_main <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank()
  )

# Mean power by regime.
p_power <- ggplot(summary_regime,
                  aes(x = as.numeric(as.character(B_tune_sec)),
                      y = mean_power,
                      group = method,
                      linetype = method,
                      shape = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_power - se_power,
                    ymax = mean_power + se_power),
                width = 5, linewidth = 0.5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_x_continuous(breaks = c(40, 80, 160, 320)) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Mean final power"
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_mean_power.png",
  plot = p_power,
  width = 10,
  height = 7,
  dpi = 300
)

# Mean valid power by regime.
p_valid_power <- ggplot(summary_regime,
                        aes(x = as.numeric(as.character(B_tune_sec)),
                            y = mean_valid_power,
                            group = method,
                            linetype = method,
                            shape = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = pmax(mean_valid_power - se_valid_power, 0),
                    ymax = mean_valid_power + se_valid_power),
                width = 5, linewidth = 0.5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_x_continuous(breaks = c(40, 80, 160, 320)) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Mean valid power"
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_mean_valid_power.png",
  plot = p_valid_power,
  width = 10,
  height = 7,
  dpi = 300
)

# Method rank heatmap.
rank_plot_df <- rank_df %>%
  mutate(
    regime = paste0("|Theta|=", theta_set_size, "\nB=", B_tune_sec),
    method = fct_rev(method)
  )

p_rank <- ggplot(rank_plot_df,
                 aes(x = B_tune_sec, y = method, fill = rank_valid_power)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = rank_valid_power), size = 5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_fill_gradient(
    low = "grey20", high = "grey90",
    limits = c(1, length(baseline_methods)),
    breaks = 1:length(baseline_methods),
    name = "Rank\n(1=best)"
  ) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = NULL
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_rank_heatmap.png",
  plot = p_rank,
  width = 10,
  height = 7,
  dpi = 300
)

# Power-runtime trade-off.
p_tradeoff <- ggplot(summary_regime,
                     aes(x = mean_c0,
                         y = mean_power,
                         shape = method,
                         linetype = method)) +
  geom_point(size = 3) +
  facet_grid(theta_set_size ~ B_tune_sec) +
  labs(
    x = expression("Mean null runtime " * E[0] * "[" * C[theta] * "]"),
    y = "Mean final power"
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_power_runtime_tradeoff.png",
  plot = p_tradeoff,
  width = 11,
  height = 8,
  dpi = 300
)

# Mean number of evaluations.
p_neval <- ggplot(summary_regime,
                  aes(x = as.numeric(as.character(B_tune_sec)),
                      y = mean_n_eval,
                      group = method,
                      linetype = method,
                      shape = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_x_continuous(breaks = c(40, 80, 160, 320)) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Mean number of evaluations"
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_n_eval.png",
  plot = p_neval,
  width = 10,
  height = 7,
  dpi = 300
)

# Power boxplots by regime.
p_box <- ggplot(df,
                aes(x = B_tune_sec, y = final_power_hat)) +
  geom_boxplot(outlier.size = 1.2, linewidth = 0.5) +
  facet_grid(theta_set_size ~ method) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Final power"
  ) +
  theme_main +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "figures/main_validonly_power_boxplots.png",
  plot = p_box,
  width = 12,
  height = 8,
  dpi = 300
)

# Valid rate plot.
p_valid_rate <- ggplot(summary_regime,
                       aes(x = as.numeric(as.character(B_tune_sec)),
                           y = mean_valid_rate,
                           group = method,
                           linetype = method,
                           shape = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_x_continuous(breaks = c(40, 80, 160, 320)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Validity rate"
  ) +
  theme_main

ggsave(
  filename = "figures/main_validonly_valid_rate.png",
  plot = p_valid_rate,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Plots saved in ./figures/\n")
