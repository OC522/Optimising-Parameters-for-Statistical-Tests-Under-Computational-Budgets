# ============================================================
# Chapter 5 methods built on algorithm3.1_runtime_fixed.R
# ============================================================
algorithm_file <- if (file.exists("algorithm3.1_runtime_fixed.R")) {
  "algorithm3.1_runtime_fixed.R"
} else {
  "~/Downloads/algorithm3.1_runtime_fixed.R"
}
source(algorithm_file)

# ------------------------------------------------------------
# Helper: final row with a separate allocated-budget diagnostic
# ------------------------------------------------------------
ch5_final_row <- function(out, method_name, seed, theta_set_id, theta_set_size,
                          B_tune_sec, allocated_budget_sec = B_tune_sec,
                          alpha_in = 0.05, N0_final = 500, N1_final = 500) {
  row <- final_row_from_tuner(out, method_name = method_name,
                              seed = seed, theta_set_id = theta_set_id,
                              theta_set_size = theta_set_size,
                              B_tune_sec = B_tune_sec,
                              alpha_in = alpha_in,
                              N0_final = N0_final, N1_final = N1_final)
  row$allocated_budget_sec <- allocated_budget_sec
  row$budget_overrun_sec <- pmax(row$tune_elapsed_sec - allocated_budget_sec, 0)
  row$budget_overrun_ratio <- ifelse(allocated_budget_sec > 0,
                                     row$tune_elapsed_sec / allocated_budget_sec,
                                     NA_real_)
  row
}

# ------------------------------------------------------------
# Ablation study for one regime / seed
# ------------------------------------------------------------
run_one_ch5_ablation <- function(seed = 1,
                                 thetas_in,
                                 theta_set_id = NA_character_,
                                 theta_set_size = length(thetas_in),
                                 B_tune_in = 160,
                                 alpha_in = 0.05,
                                 N0_final = 500,
                                 N1_final = 500,
                                 lambda_sh = 0.50,
                                 r_eval = 0.25,
                                 sh_eta = 2,
                                 sh_lcb_conf = 0.90,
                                 sa_T0 = 0.20,
                                 sa_cooling = 0.95) {
  set.seed(seed)
  B_sh <- lambda_sh * B_tune_in
  B_sa <- (1 - lambda_sh) * B_tune_in

  out_sh_full <- successive_halving_timebudgeted(
    thetas = thetas_in, B_tune_sec = B_tune_in,
    eta = sh_eta, alpha = alpha_in, rH1 = rH1_shift005,
    lcb_conf = sh_lcb_conf, seed = seed)

  out_sa_full <- tune_sa_valid_graph(
    thetas = thetas_in, B_tune_sec = B_tune_in,
    r_eval = r_eval, T0 = sa_T0, cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed,
    method_label = "SA_full")

  out_hybrid <- tune_hybrid_sh_sa(
    thetas = thetas_in, B_tune_sec = B_tune_in,
    lambda_sh = lambda_sh, sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
    r_eval = r_eval, sa_T0 = sa_T0, sa_cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed,
    restrict_to_sh_survivors = FALSE)

  out_sh_half <- successive_halving_timebudgeted(
    thetas = thetas_in, B_tune_sec = B_sh,
    eta = sh_eta, alpha = alpha_in, rH1 = rH1_shift005,
    lcb_conf = sh_lcb_conf, seed = seed)

  out_sa_random_half <- tune_sa_valid_graph(
    thetas = thetas_in, B_tune_sec = B_sa,
    r_eval = r_eval, T0 = sa_T0, cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed + 10000L,
    method_label = "SA_random_half")

  # SH initialisation is paid for, but this row reports the local-refinement
  # result alone.  This distinguishes it from Hybrid, which returns the best
  # incumbent over both stages.
  init_obj <- choose_hybrid_init_from_sh(out_sh_half, thetas_in)
  B_sa_after_init <- max(0, B_tune_in - out_sh_half$tune_elapsed)
  out_sa_from_sh_only <- tune_sa_valid_graph(
    thetas = thetas_in, B_tune_sec = B_sa_after_init,
    init_idx = init_obj$init_idx,
    r_eval = r_eval, T0 = sa_T0, cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed,
    method_label = "SA_from_SH_half")
  out_sa_from_sh_half <- out_sa_from_sh_only
  out_sa_from_sh_half$method <- "SA_from_SH_half"
  out_sa_from_sh_half$tune_elapsed <- out_sh_half$tune_elapsed + out_sa_from_sh_only$tune_elapsed
  out_sa_from_sh_half$n_theta_eval <- out_sh_half$n_theta_eval + out_sa_from_sh_only$n_theta_eval
  out_sa_from_sh_half$n_timed_out <- out_sh_half$n_timed_out + out_sa_from_sh_only$n_timed_out

  out_hybrid_restricted <- tune_hybrid_sh_sa(
    thetas = thetas_in, B_tune_sec = B_tune_in,
    lambda_sh = lambda_sh, sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
    r_eval = r_eval, sa_T0 = sa_T0, sa_cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed,
    restrict_to_sh_survivors = TRUE)

  rows <- list(
    ch5_final_row(out_sh_full, "SH_full", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_tune_in, alpha_in, N0_final, N1_final),
    ch5_final_row(out_sa_full, "SA_full", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_tune_in, alpha_in, N0_final, N1_final),
    ch5_final_row(out_hybrid, "Hybrid", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_tune_in, alpha_in, N0_final, N1_final),
    ch5_final_row(out_sh_half, "SH_half", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_sh, alpha_in, N0_final, N1_final),
    ch5_final_row(out_sa_random_half, "SA_random_half", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_sa, alpha_in, N0_final, N1_final),
    ch5_final_row(out_sa_from_sh_half, "SA_from_SH_half", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_tune_in, alpha_in, N0_final, N1_final),
    ch5_final_row(out_hybrid_restricted, "Hybrid_restricted", seed, theta_set_id, theta_set_size,
                  B_tune_in, B_tune_in, alpha_in, N0_final, N1_final)
  )
  dplyr::bind_rows(rows)
}

run_ch5_ablation_grid <- function(all_thetas = full_thetas,
                                  candidate_sizes = c(128, 512),
                                  budget_grid = c(80, 160, 320),
                                  seeds = 1:3,
                                  alpha_in = 0.05,
                                  N0_final = 500,
                                  N1_final = 500,
                                  lambda_sh = 0.50,
                                  r_eval = 0.25,
                                  sh_eta = 2,
                                  sh_lcb_conf = 0.90,
                                  sa_T0 = 0.20,
                                  sa_cooling = 0.95,
                                  output_csv = "results_ch5_hybrid_ablation_runtime_fixed.csv") {
  rows <- list(); id <- 0L
  for (M in candidate_sizes) {
    for (s in seeds) {
      subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
      theta_set_id <- paste0("CH5_M", M, "_seed", s)
      for (B_tune in budget_grid) {
        cat("CH5 ablation: |Theta| =", M, ", B =", B_tune,
            ", seed =", s, "\n")
        id <- id + 1L
        rows[[id]] <- run_one_ch5_ablation(
          seed = s, thetas_in = subset_obj$thetas,
          theta_set_id = theta_set_id, theta_set_size = M,
          B_tune_in = B_tune, alpha_in = alpha_in,
          N0_final = N0_final, N1_final = N1_final,
          lambda_sh = lambda_sh, r_eval = r_eval,
          sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
          sa_T0 = sa_T0, sa_cooling = sa_cooling)
      }
    }
  }
  out_df <- dplyr::bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

# ------------------------------------------------------------
# Lambda sensitivity
# ------------------------------------------------------------
run_ch5_lambda_sensitivity <- function(all_thetas = full_thetas,
                                       candidate_sizes = c(128, 512),
                                       budget_grid = c(80, 160, 320),
                                       seeds = 1:3,
                                       lambda_grid = c(0.30, 0.50, 0.70),
                                       alpha_in = 0.05,
                                       N0_final = 500,
                                       N1_final = 500,
                                       r_eval = 0.25,
                                       sh_eta = 2,
                                       sh_lcb_conf = 0.90,
                                       sa_T0 = 0.20,
                                       sa_cooling = 0.95,
                                       output_csv = "results_ch5_lambda_sensitivity_runtime_fixed.csv") {
  rows <- list(); id <- 0L
  for (M in candidate_sizes) {
    for (s in seeds) {
      subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
      theta_set_id <- paste0("CH5LAMBDA_M", M, "_seed", s)
      for (B_tune in budget_grid) {
        for (lambda_val in lambda_grid) {
          cat("CH5 lambda: |Theta| =", M, ", B =", B_tune,
              ", seed =", s, ", lambda =", lambda_val, "\n")
          out <- tune_hybrid_sh_sa(
            thetas = subset_obj$thetas, B_tune_sec = B_tune,
            lambda_sh = lambda_val, sh_eta = sh_eta,
            sh_lcb_conf = sh_lcb_conf, r_eval = r_eval,
            sa_T0 = sa_T0, sa_cooling = sa_cooling,
            alpha = alpha_in, rH1 = rH1_shift005, seed = s)
          id <- id + 1L
          rows[[id]] <- ch5_final_row(out, paste0("Hybrid_lambda_", lambda_val),
                                      s, theta_set_id, M, B_tune, B_tune,
                                      alpha_in, N0_final, N1_final)
        }
      }
    }
  }
  out_df <- dplyr::bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

# ------------------------------------------------------------
# Wrapper to run all Chapter 5 experiments
# ------------------------------------------------------------
run_ch5_all_runtime_fixed <- function() {
  main <- run_experiment_grid_with_hybrid(
    all_thetas = full_thetas,
    candidate_sizes = c(64, 128, 256, 512),
    budget_grid = c(40, 80, 160, 320),
    seeds = 1:5,
    N0_final = 500,
    N1_final = 500,
    output_csv = "results_main_with_hybrid_runtime_fixed.csv")

  abl <- run_ch5_ablation_grid(
    all_thetas = full_thetas,
    candidate_sizes = c(128),
    budget_grid = c(80, 160, 320),
    seeds = 1:3,
    N0_final = 500,
    N1_final = 500,
    output_csv = "results_ch5_hybrid_ablation_runtime_fixed.csv")

  lam <- run_ch5_lambda_sensitivity(
    all_thetas = full_thetas,
    candidate_sizes = c(128, 512),
    budget_grid = c(80, 160, 320),
    seeds = 1:3,
    N0_final = 500,
    N1_final = 500,
    output_csv = "results_ch5_lambda_sensitivity_runtime_fixed.csv")

  list(main = main, ablation = abl, lambda = lam)
}
