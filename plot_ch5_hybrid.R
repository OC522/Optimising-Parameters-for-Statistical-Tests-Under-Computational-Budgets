# Older Chapter 5 plotting script.

library(dplyr)
library(ggplot2)
library(forcats)

if (!dir.exists("figures")) dir.create("figures")

theme_main <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

# Shapes used across the method sets.
shape_values_5 <- c(
  "Uniform" = 16,
  "Random"  = 17,
  "SA"      = 15,
  "SH"      = 3,
  "Hybrid"  = 7
)

shape_values_7 <- c(
  "SH_full"           = 16,
  "SA_full"           = 17,
  "Hybrid"            = 15,
  "SH_half"           = 3,
  "SA_random_half"    = 7,
  "SA_from_SH_half"   = 8,
  "Hybrid_restricted" = 4
)

# Standard error with a small NA guard.
safe_se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(0)
  sd(x) / sqrt(length(x))
}

# Main comparison with the hybrid.
main_file <- if (file.exists("results_main_with_hybrid_fair.csv")) {
  "results_main_with_hybrid_fair.csv"
} else {
  "results_main_validonly.csv"
}

main_df <- read.csv(main_file, stringsAsFactors = FALSE)

main_methods <- c("Uniform", "Random", "SA", "SH", "Hybrid")

main_df <- main_df %>%
  filter(method %in% main_methods) %>%
  mutate(
    method = factor(method, levels = main_methods),
    theta_set_size = factor(theta_set_size, levels = c(64, 128, 256, 512)),
    B_tune_sec = factor(B_tune_sec, levels = c(40, 80, 160, 320))
  )

main_sum <- main_df %>%
  group_by(theta_set_size, B_tune_sec, method) %>%
  summarise(
    n_valid_power = sum(!is.na(valid_power)),
    mean_valid_power = mean(valid_power, na.rm = TRUE),
    se_valid_power = safe_se(valid_power),
    mean_power = mean(final_power_hat, na.rm = TRUE),
    mean_type1 = mean(final_type1_hat, na.rm = TRUE),
    mean_valid_rate = mean(valid_flag, na.rm = TRUE),
    mean_c0 = mean(final_c0_mean, na.rm = TRUE),
    mean_n_eval = mean(n_theta_eval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    n_valid_power > 0,
    is.finite(mean_valid_power),
    is.finite(se_valid_power)
  )

p_main <- ggplot(
  main_sum,
  aes(
    x = as.numeric(as.character(B_tune_sec)),
    y = mean_valid_power,
    group = method,
    linetype = method,
    shape = method
  )
) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_valid_power - se_valid_power, 0),
      ymax = pmin(mean_valid_power + se_valid_power, 1)
    ),
    width = 5,
    linewidth = 0.5,
    na.rm = TRUE
  ) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_x_continuous(breaks = c(40, 80, 160, 320)) +
  scale_shape_manual(values = shape_values_5, drop = FALSE) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = "Mean valid power"
  ) +
  theme_main

ggsave(
  "figures/ch5_main_mean_valid_power_with_hybrid.png",
  p_main,
  width = 10,
  height = 7,
  dpi = 300
)

rank_df <- main_sum %>%
  group_by(theta_set_size, B_tune_sec) %>%
  arrange(desc(mean_valid_power), .by_group = TRUE) %>%
  mutate(rank_valid_power = row_number()) %>%
  ungroup() %>%
  mutate(method = fct_rev(method))

p_rank <- ggplot(rank_df, aes(x = B_tune_sec, y = method, fill = rank_valid_power)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = rank_valid_power), size = 5) +
  facet_wrap(~ theta_set_size, nrow = 2) +
  scale_fill_gradient(
    low = "grey20",
    high = "grey90",
    breaks = 1:length(main_methods),
    name = "Rank\n(1=best)"
  ) +
  labs(
    x = expression(B[tune] ~ "(seconds)"),
    y = NULL
  ) +
  theme_main

ggsave(
  "figures/ch5_main_rank_heatmap_with_hybrid.png",
  p_rank,
  width = 10,
  height = 7,
  dpi = 300
)

# Ablation plot.
if (file.exists("results_ch5_hybrid_ablation.csv")) {
  abl <- read.csv("results_ch5_hybrid_ablation.csv", stringsAsFactors = FALSE)
  
  abl_methods <- c(
    "SH_full",
    "SA_full",
    "Hybrid",
    "SH_half",
    "SA_random_half",
    "SA_from_SH_half",
    "Hybrid_restricted"
  )
  
  abl <- abl %>%
    filter(method %in% abl_methods) %>%
    mutate(
      method = factor(method, levels = abl_methods),
      theta_set_size = factor(theta_set_size),
      B_tune_sec = factor(B_tune_sec)
    )
  
  abl_sum <- abl %>%
    group_by(theta_set_size, B_tune_sec, method) %>%
    summarise(
      n_valid_power = sum(!is.na(valid_power)),
      mean_valid_power = mean(valid_power, na.rm = TRUE),
      se_valid_power = safe_se(valid_power),
      mean_valid_rate = mean(valid_flag, na.rm = TRUE),
      mean_power = mean(final_power_hat, na.rm = TRUE),
      mean_c0 = mean(final_c0_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      n_valid_power > 0,
      is.finite(mean_valid_power)
    )
  
  p_abl <- ggplot(
    abl_sum,
    aes(
      x = B_tune_sec,
      y = mean_valid_power,
      group = method,
      linetype = method,
      shape = method
    )
  ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.4, na.rm = TRUE) +
    facet_wrap(~ theta_set_size, nrow = 1) +
    scale_shape_manual(values = shape_values_7, drop = FALSE) +
    labs(
      x = expression(B[tune] ~ "(seconds)"),
      y = "Mean valid power"
    ) +
    theme_main
  
  ggsave(
    "figures/ch5_hybrid_ablation_valid_power.png",
    p_abl,
    width = 11,
    height = 6,
    dpi = 300
  )
}

# Lambda sensitivity plot.
if (file.exists("results_ch5_lambda_sensitivity.csv")) {
  lam <- read.csv("results_ch5_lambda_sensitivity.csv", stringsAsFactors = FALSE)
  
  lam <- lam %>%
    filter(grepl("^Hybrid_lambda_", method)) %>%
    mutate(
      lambda_sh = as.numeric(gsub("Hybrid_lambda_", "", method)),
      theta_set_size = factor(theta_set_size),
      B_tune_sec = factor(B_tune_sec)
    ) %>%
    filter(!is.na(lambda_sh))
  
  lam_sum <- lam %>%
    group_by(theta_set_size, B_tune_sec, lambda_sh) %>%
    summarise(
      n_valid_power = sum(!is.na(valid_power)),
      mean_valid_power = mean(valid_power, na.rm = TRUE),
      mean_valid_rate = mean(valid_flag, na.rm = TRUE),
      mean_power = mean(final_power_hat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      n_valid_power > 0,
      is.finite(mean_valid_power)
    )
  
  p_lam <- ggplot(
    lam_sum,
    aes(
      x = factor(lambda_sh),
      y = mean_valid_power,
      group = B_tune_sec,
      linetype = B_tune_sec,
      shape = B_tune_sec
    )
  ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.5, na.rm = TRUE) +
    facet_wrap(~ theta_set_size, nrow = 1) +
    labs(
      x = expression(lambda),
      y = "Mean valid power",
      linetype = expression(B[tune]),
      shape = expression(B[tune])
    ) +
    theme_main
  
  ggsave(
    "figures/ch5_lambda_sensitivity.png",
    p_lam,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# Console summaries.
message("Chapter 5 figures written to ./figures/")
message("Generated files:")
print(list.files("figures", pattern = "^ch5_.*\\.png$", full.names = TRUE))
