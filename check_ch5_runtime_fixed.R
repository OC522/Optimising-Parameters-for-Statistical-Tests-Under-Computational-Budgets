# ============================================================
# Diagnostics for runtime-fixed Chapter 5 results
# ============================================================
algorithm_file <- if (file.exists("algorithm3.1_runtime_fixed.R")) {
  "algorithm3.1_runtime_fixed.R"
} else {
  "~/Downloads/algorithm3.1_runtime_fixed.R"
}
source(algorithm_file)

required_files <- c(
  main = "results_main_with_hybrid_runtime_fixed.csv",
  ablation = "results_ch5_hybrid_ablation_runtime_fixed.csv",
  lambda = "results_ch5_lambda_sensitivity_runtime_fixed.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing runtime-fixed result file(s): ",
       paste(missing_files, collapse = ", "), call. = FALSE)
}

main <- read.csv(required_files["main"], stringsAsFactors = FALSE)
abl  <- read.csv(required_files["ablation"], stringsAsFactors = FALSE)
lam  <- read.csv(required_files["lambda"], stringsAsFactors = FALSE)

cat("\nResult-file dimensions\n")
cat("MAIN rows:", nrow(main), "(expected 400)\n")
cat("ABLATION rows:", nrow(abl), "(expected 63 for |Theta|=128 only)\n")
cat("LAMBDA rows:", nrow(lam), "(expected 54)\n")

if (!is.null(main)) {
  cat("\nMAIN budget diagnostics\n")
  print(budget_diagnostics(main))
}
if (!is.null(abl)) {
  cat("\nABLATION budget diagnostics\n")
  print(budget_diagnostics(abl))
}
if (!is.null(lam)) {
  cat("\nLAMBDA budget diagnostics\n")
  print(budget_diagnostics(lam))
}

summarise_methods <- function(df) {
  df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(mean_valid_power = mean(valid_power, na.rm = TRUE),
                     mean_power = mean(final_power_hat, na.rm = TRUE),
                     mean_type1 = mean(final_type1_hat, na.rm = TRUE),
                     mean_valid_rate = mean(valid_flag, na.rm = TRUE),
                     mean_c0 = mean(final_c0_mean, na.rm = TRUE),
                     mean_n_eval = mean(n_theta_eval, na.rm = TRUE),
                     mean_n_timed_out = mean(n_timed_out, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_valid_power))
}

if (!is.null(main)) {
  cat("\nMAIN method summary\n")
  print(summarise_methods(main))
}

cat("\nMAIN thesis-number checks\n")
main_sum <- main %>%
  dplyr::group_by(theta_set_size, B_tune_sec, method) %>%
  dplyr::summarise(mean_valid_power = mean(valid_power, na.rm = TRUE),
                   .groups = "drop")
rank_df <- main_sum %>%
  dplyr::group_by(theta_set_size, B_tune_sec) %>%
  dplyr::arrange(dplyr::desc(mean_valid_power), .by_group = TRUE) %>%
  dplyr::mutate(rank_valid_power = dplyr::row_number()) %>%
  dplyr::ungroup()
print(rank_df %>%
        dplyr::filter(method == "Hybrid") %>%
        dplyr::arrange(theta_set_size, B_tune_sec))
cat("Hybrid rank-1 regimes:",
    sum(rank_df$method == "Hybrid" & rank_df$rank_valid_power == 1), "\n")
cat("Hybrid top-2 regimes:",
    sum(rank_df$method == "Hybrid" & rank_df$rank_valid_power <= 2), "\n")

main_wide <- reshape(as.data.frame(main_sum),
                     idvar = c("theta_set_size", "B_tune_sec"),
                     timevar = "method", direction = "wide")
cat("Hybrid-SA mean difference:",
    mean(main_wide$mean_valid_power.Hybrid - main_wide$mean_valid_power.SA), "\n")
cat("Hybrid-SH mean difference:",
    mean(main_wide$mean_valid_power.Hybrid - main_wide$mean_valid_power.SH), "\n")
cat("Hybrid-Uniform mean difference:",
    mean(main_wide$mean_valid_power.Hybrid - main_wide$mean_valid_power.Uniform), "\n")
cat("Hybrid-Random mean difference:",
    mean(main_wide$mean_valid_power.Hybrid - main_wide$mean_valid_power.Random), "\n")

cat("\nABLATION table for |Theta|=128\n")
print(abl %>%
        dplyr::filter(theta_set_size == 128) %>%
        dplyr::group_by(method, B_tune_sec) %>%
        dplyr::summarise(mean_valid_power = round(mean(valid_power, na.rm = TRUE), 3),
                         .groups = "drop") %>%
        tidyr::pivot_wider(names_from = B_tune_sec, values_from = mean_valid_power) %>%
        dplyr::arrange(match(method, c("SH_full", "SA_full", "Hybrid",
                                       "SH_half", "SA_random_half",
                                       "SA_from_SH_half",
                                       "Hybrid_restricted"))))

cat("\nLAMBDA thesis-number checks\n")
lam2 <- lam %>%
  dplyr::mutate(lambda_sh = as.numeric(gsub("Hybrid_lambda_", "", method)))
print(lam2 %>%
        dplyr::group_by(lambda_sh) %>%
        dplyr::summarise(mean_valid_power = mean(valid_power, na.rm = TRUE),
                         valid_rate = mean(valid_flag, na.rm = TRUE),
                         .groups = "drop"))
lam_rank <- lam2 %>%
  dplyr::group_by(theta_set_size, B_tune_sec, lambda_sh) %>%
  dplyr::summarise(mean_valid_power = mean(valid_power, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::group_by(theta_set_size, B_tune_sec) %>%
  dplyr::arrange(dplyr::desc(mean_valid_power), .by_group = TRUE) %>%
  dplyr::mutate(rank_valid_power = dplyr::row_number()) %>%
  dplyr::ungroup()
print(lam_rank %>%
        dplyr::filter(rank_valid_power == 1) %>%
        dplyr::count(lambda_sh, name = "rank1_regimes"))
