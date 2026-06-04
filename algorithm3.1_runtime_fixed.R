# =========================================
# Sequential KS tuning under runtime budgets
# Runtime-fixed version for Chapters 3--5
#
# Key fixes relative to the earlier algorithm3.1.R:
#   1. seq_ks_once() is lazy stage-by-stage, not full n[K] upfront.
#   2. tuning-time evaluations receive an explicit deadline and sample guard.
#   3. min_runs is a target, not permission to exceed the runtime budget.
#   4. SA moves on the fixed candidate set via a kNN graph, preserving fair search spaces.
#   5. SH uses complete-round logic and budget-aware per-candidate evaluation.
#   6. final reporting includes budget diagnostics and deterministic final seeds.
# =========================================

set.seed(1)
alpha <- 0.05

suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------------------
# Alternative used under H1
# -----------------------------------------
# Function name retained for compatibility with earlier drafts; the current
# benchmark alternative is N(0.03, 1), matching the thesis text.
rH1_shift005 <- function(n) rnorm(n, mean = 0.03, sd = 1)
# Alias retained for compatibility with older Chapter 5 files.
rH1_shift003 <- rH1_shift005

# -----------------------------------------
# Basic utilities
# -----------------------------------------
now_elapsed <- function() proc.time()[["elapsed"]]

remaining_budget <- function(t0, B) {
  max(0, as.numeric(B) - (now_elapsed() - t0))
}

alpha_spend_geom <- function(alpha, K, rho) {
  w <- rho^(0:(K - 1))
  alpha * w / sum(w)
}

stage_sizes <- function(n0, gamma, K) {
  as.integer(ceiling(n0 * gamma^(0:(K - 1))))
}

wilson_lower <- function(x, n, conf = 0.90) {
  if (!is.finite(x) || !is.finite(n) || n <= 0) return(0)
  z <- qnorm(conf)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- (z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))) / denom
  max(0, center - half)
}

wilson_ci <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- (z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))) / denom
  c(max(0, center - half), min(1, center + half))
}

theta_id_string <- function(theta) {
  paste(theta$n0, theta$gamma, theta$K, theta$rho, theta$u, sep = "|")
}

stable_seed <- function(...) {
  s <- paste(..., sep = "|")
  z <- utf8ToInt(s)
  v <- sum(z * seq_along(z)) %% 2147483000L
  as.integer(max(1L, v))
}

# -----------------------------------------
# Fast asymptotic one-sample KS p-value
# -----------------------------------------
kolmogorov_sf <- function(z, tol = 1e-12, max_terms = 100L) {
  if (!is.finite(z)) return(NA_real_)
  if (z <= 0) return(1)
  j <- seq_len(max_terms)
  terms <- 2 * (-1)^(j - 1) * exp(-2 * (j^2) * (z^2))
  # Truncate when terms are negligible.
  small <- which(abs(terms) < tol)
  if (length(small) > 0) terms <- terms[seq_len(small[1])]
  p <- sum(terms)
  min(1, max(0, p))
}

fast_ks_pvalue <- function(x) {
  n <- length(x)
  if (n <= 0) return(NA_real_)
  xs <- sort(x)
  Fx <- pnorm(xs)
  i <- seq_len(n)
  D_plus <- max(i / n - Fx)
  D_minus <- max(Fx - (i - 1) / n)
  D <- max(D_plus, D_minus)
  z <- (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D
  kolmogorov_sf(z)
}

# Conservative sample guard used only during tuning.  The goal is not to
# approximate runtime perfectly, but to prevent a candidate from entering a
# million-sample stage when the per-candidate budget is a fraction of a second.
max_n_for_tuning_budget <- function(seconds,
                                    min_n = 2000L,
                                    base_n = 15000L,
                                    per_sec = 60000L,
                                    hard_max = 120000L) {
  if (!is.finite(seconds) || seconds <= 0) return(0L)
  as.integer(max(min_n, min(hard_max, base_n + per_sec * seconds)))
}

# -----------------------------------------
# One run of sequential KS
# -----------------------------------------
seq_ks_once <- function(theta,
                        alpha = 0.05,
                        scenario = c("H0", "H1"),
                        rH1 = rH1_shift005,
                        deadline = Inf,
                        max_n_total = Inf,
                        use_fast_ks = TRUE) {
  scenario <- match.arg(scenario)

  n0    <- theta$n0
  gamma <- theta$gamma
  K     <- theta$K
  rho   <- theta$rho
  u     <- theta$u

  a <- alpha_spend_geom(alpha, K, rho)
  n <- stage_sizes(n0, gamma, K)

  t0 <- now_elapsed()
  x <- numeric(0)
  n_prev <- 0L
  p <- NA_real_

  for (k in seq_len(K)) {
    if (now_elapsed() >= deadline) {
      return(list(reject = NA_integer_, stage = k, p = p,
                  cpu = now_elapsed() - t0, aborted = TRUE,
                  abort_reason = "deadline"))
    }
    if (n[k] > max_n_total) {
      return(list(reject = NA_integer_, stage = k, p = p,
                  cpu = now_elapsed() - t0, aborted = TRUE,
                  abort_reason = "sample_guard"))
    }

    need <- as.integer(n[k] - n_prev)
    if (need > 0L) {
      x_new <- if (scenario == "H0") rnorm(need) else rH1(need)
      x <- c(x, x_new)
      n_prev <- n[k]
    }

    p <- if (use_fast_ks) {
      fast_ks_pvalue(x)
    } else {
      suppressWarnings(stats::ks.test(x, "pnorm", exact = FALSE)$p.value)
    }

    if (now_elapsed() >= deadline) {
      return(list(reject = NA_integer_, stage = k, p = p,
                  cpu = now_elapsed() - t0, aborted = TRUE,
                  abort_reason = "deadline"))
    }

    if (is.finite(p) && p <= a[k]) {
      return(list(reject = 1L, stage = k, p = p,
                  cpu = now_elapsed() - t0, aborted = FALSE,
                  abort_reason = NA_character_))
    }
    if (is.finite(p) && p >= u) {
      return(list(reject = 0L, stage = k, p = p,
                  cpu = now_elapsed() - t0, aborted = FALSE,
                  abort_reason = NA_character_))
    }
  }

  list(reject = 0L, stage = K, p = p,
       cpu = now_elapsed() - t0, aborted = FALSE,
       abort_reason = NA_character_)
}

# -----------------------------------------
# Budgeted Monte Carlo evaluation of one scenario
# -----------------------------------------
eval_theta_timebudget <- function(theta,
                                  time_budget_sec,
                                  alpha = 0.05,
                                  scenario = c("H1", "H0"),
                                  rH1 = rH1_shift005,
                                  min_runs = 1L,
                                  deadline = NULL,
                                  max_n_total = NULL,
                                  use_fast_ks = TRUE) {
  scenario <- match.arg(scenario)
  min_runs <- as.integer(min_runs)
  start <- now_elapsed()
  if (is.null(deadline)) deadline <- start + max(0, time_budget_sec)
  if (is.null(max_n_total)) max_n_total <- max_n_for_tuning_budget(time_budget_sec)

  n_run <- 0L
  n_rej <- 0L
  cpu_sum <- 0
  n_aborted <- 0L
  abort_reasons <- character(0)

  repeat {
    if ((now_elapsed() - start) >= time_budget_sec) break
    if (now_elapsed() >= deadline) break

    out <- seq_ks_once(theta, alpha = alpha, scenario = scenario, rH1 = rH1,
                       deadline = deadline, max_n_total = max_n_total,
                       use_fast_ks = use_fast_ks)

    if (isTRUE(out$aborted)) {
      n_aborted <- n_aborted + 1L
      abort_reasons <- c(abort_reasons, out$abort_reason)
      break
    }

    n_run <- n_run + 1L
    n_rej <- n_rej + as.integer(out$reject)
    cpu_sum <- cpu_sum + out$cpu
  }

  enough_runs <- n_run >= min_runs
  rate <- if (n_run > 0L) n_rej / n_run else NA_real_

  list(
    rate = if (enough_runs) rate else NA_real_,
    n_run = n_run,
    n_rej = n_rej,
    avg_cpu_per_run = if (n_run > 0L) cpu_sum / n_run else NA_real_,
    elapsed = now_elapsed() - start,
    enough_runs = enough_runs,
    timed_out = !enough_runs || n_aborted > 0L || now_elapsed() >= deadline,
    n_aborted = n_aborted,
    abort_reasons = paste(unique(abort_reasons), collapse = ",")
  )
}

# -----------------------------------------
# Candidate set Theta (full grid)
# -----------------------------------------
make_candidates <- function() {
  grid <- expand.grid(
    n0    = c(40, 60, 80, 100, 120, 160, 200, 260, 320),
    gamma = c(1.10, 1.18, 1.25, 1.35, 1.50, 1.65, 1.85, 2.10, 2.40),
    K     = c(5, 6, 7, 8),
    rho   = c(0.35, 0.45, 0.60, 0.75, 0.90, 0.97),
    u     = c(0.80, 0.86, 0.90, 0.92, 0.95, 0.98),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    grid = grid,
    thetas = lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, , drop = FALSE])),
    levels = list(
      n0 = sort(unique(grid$n0)),
      gamma = sort(unique(grid$gamma)),
      K = sort(unique(grid$K)),
      rho = sort(unique(grid$rho)),
      u = sort(unique(grid$u))
    )
  )
}

cand <- make_candidates()
full_grid_df <- cand$grid
full_thetas <- cand$thetas
full_theta_levels <- cand$levels
cat("Full candidate count |Theta| =", length(full_thetas), "\n")

# -----------------------------------------
# Helpers for subsets / indexing / graph neighbours
# -----------------------------------------
make_levels_from_thetas <- function(thetas) {
  df <- do.call(rbind, lapply(thetas, function(th) {
    data.frame(n0 = th$n0, gamma = th$gamma, K = th$K,
               rho = th$rho, u = th$u)
  }))
  list(n0 = sort(unique(df$n0)), gamma = sort(unique(df$gamma)),
       K = sort(unique(df$K)), rho = sort(unique(df$rho)),
       u = sort(unique(df$u)))
}

sample_candidate_subset <- function(all_thetas, M, seed) {
  stopifnot(M <= length(all_thetas))
  set.seed(seed)
  idx <- sample.int(length(all_thetas), M, replace = FALSE)
  thetas_sub <- all_thetas[idx]
  list(idx = idx, thetas = thetas_sub, levels = make_levels_from_thetas(thetas_sub))
}

theta_to_index <- function(theta, levels) {
  c(n0 = match(theta$n0, levels$n0),
    gamma = match(theta$gamma, levels$gamma),
    K = match(theta$K, levels$K),
    rho = match(theta$rho, levels$rho),
    u = match(theta$u, levels$u))
}

index_to_theta <- function(idx_vec, levels) {
  list(n0 = levels$n0[idx_vec["n0"]],
       gamma = levels$gamma[idx_vec["gamma"]],
       K = levels$K[idx_vec["K"]],
       rho = levels$rho[idx_vec["rho"]],
       u = levels$u[idx_vec["u"]])
}

build_theta_lookup <- function(thetas) {
  ids <- vapply(thetas, theta_id_string, character(1))
  structure(seq_along(thetas), names = ids)
}

propose_neighbor <- function(idx_vec, levels) {
  dims <- names(idx_vec)
  repeat {
    d <- sample(dims, 1)
    max_idx <- c(n0 = length(levels$n0), gamma = length(levels$gamma),
                 K = length(levels$K), rho = length(levels$rho), u = length(levels$u))
    step <- sample(c(-1L, 1L), 1)
    proposed <- idx_vec[d] + step
    if (proposed >= 1L && proposed <= max_idx[d]) {
      idx_new <- idx_vec
      idx_new[d] <- proposed
      return(idx_new)
    }
  }
}

thetas_to_df <- function(thetas) {
  do.call(rbind, lapply(thetas, function(th) {
    data.frame(n0 = th$n0, gamma = th$gamma, K = th$K,
               rho = th$rho, u = th$u)
  }))
}

encode_thetas_scaled <- function(thetas, levels = make_levels_from_thetas(thetas)) {
  df <- thetas_to_df(thetas)
  enc <- cbind(
    n0 = match(df$n0, levels$n0),
    gamma = match(df$gamma, levels$gamma),
    K = match(df$K, levels$K),
    rho = match(df$rho, levels$rho),
    u = match(df$u, levels$u)
  )
  enc <- apply(enc, 2, function(v) {
    rng <- max(v) - min(v)
    if (rng == 0) rep(0, length(v)) else (v - min(v)) / rng
  })
  as.matrix(enc)
}

build_knn_neighbors <- function(thetas, levels = make_levels_from_thetas(thetas), k_neighbors = 8L) {
  M <- length(thetas)
  if (M <= 1L) return(vector("list", M))
  X <- encode_thetas_scaled(thetas, levels)
  k <- min(as.integer(k_neighbors), M - 1L)
  out <- vector("list", M)
  for (i in seq_len(M)) {
    d <- rowSums((X - matrix(X[i, ], nrow = M, ncol = ncol(X), byrow = TRUE))^2)
    out[[i]] <- order(d)[2:(k + 1L)]
  }
  out
}

propose_neighbor_graph <- function(cur_idx, neighbors) {
  nb <- neighbors[[cur_idx]]
  if (length(nb) == 0L) return(cur_idx)
  sample(nb, 1L)
}

last_survivors_from_sh <- function(out_sh) {
  if (is.null(out_sh$history) || length(out_sh$history) == 0L) return(out_sh$best_idx)
  h <- out_sh$history[[length(out_sh$history)]]
  if (!is.null(h$survivors) && length(h$survivors) > 0L) return(h$survivors)
  out_sh$best_idx
}

choose_hybrid_init_from_sh <- function(out_sh, thetas) {
  init_idx <- out_sh$best_idx
  if (!is.null(out_sh$history) && length(out_sh$history) > 0) {
    last_hist <- out_sh$history[[length(out_sh$history)]]
    tab <- last_hist$table
    if (!is.null(tab) && nrow(tab) > 0) {
      rank_col <- if ("rank_loss" %in% names(tab)) "rank_loss" else "loss"
      valid_tab <- tab[is.finite(tab[[rank_col]]), , drop = FALSE]
      if (nrow(valid_tab) > 0) {
        valid_tab <- valid_tab[order(valid_tab[[rank_col]], -valid_tab$power_n), , drop = FALSE]
        init_idx <- valid_tab$idx[1]
      }
    }
  }
  if (!is.finite(init_idx) || init_idx < 1L) init_idx <- 1L
  list(init_idx = init_idx, init_theta = thetas[[init_idx]])
}

# -----------------------------------------
# Final evaluation
# -----------------------------------------
final_eval_fixedN <- function(theta, N, alpha = 0.05,
                              scenario = c("H0", "H1"),
                              rH1 = rH1_shift005,
                              use_fast_ks = TRUE) {
  scenario <- match.arg(scenario)
  rej <- integer(N)
  cpu <- numeric(N)
  stage <- integer(N)
  aborted <- logical(N)

  for (i in seq_len(N)) {
    out1 <- seq_ks_once(theta, alpha = alpha, scenario = scenario, rH1 = rH1,
                        deadline = Inf, max_n_total = Inf,
                        use_fast_ks = use_fast_ks)
    rej[i] <- if (isTRUE(out1$aborted)) 0L else as.integer(out1$reject)
    cpu[i] <- out1$cpu
    stage[i] <- out1$stage
    aborted[i] <- isTRUE(out1$aborted)
  }

  rate <- mean(rej)
  ci <- wilson_ci(sum(rej), N, conf = 0.95)
  list(N = N, rate = rate, ci_low = ci[1], ci_high = ci[2],
       cpu_mean = mean(cpu), cpu_median = median(cpu),
       cpu_q90 = as.numeric(quantile(cpu, 0.9)),
       cpu_q95 = as.numeric(quantile(cpu, 0.95)),
       stage_mean = mean(stage), stage_median = median(stage),
       n_aborted = sum(aborted),
       raw = list(rej = rej, cpu = cpu, stage = stage, aborted = aborted))
}

# -----------------------------------------
# Common validity-aware oracle
# -----------------------------------------
evaluate_theta_valid <- function(theta, r,
                                 alpha = 0.05,
                                 rH1 = rH1_shift005,
                                 min_h0_runs = 2L,
                                 min_h1_runs = 6L,
                                 hard_timeout = NULL,
                                 use_fast_ks = TRUE) {
  min_h0_runs <- as.integer(min_h0_runs)
  min_h1_runs <- as.integer(min_h1_runs)

  if (is.null(hard_timeout)) hard_timeout <- min(2.5, 4.0 * r + 0.10)
  hard_timeout <- as.numeric(hard_timeout)

  if (!is.finite(hard_timeout) || hard_timeout <= 0.05) {
    return(list(valid = FALSE, type1_hat = NA_real_, type1_n = 0L,
                cpu0_hat = NA_real_, power_hat = NA_real_, power_n = 0L,
                power_rej = 0L, loss = Inf, timed_out = TRUE,
                h0_timed_out = TRUE, h1_timed_out = TRUE))
  }

  t_eval <- now_elapsed()
  deadline <- t_eval + hard_timeout
  max_n_total <- max_n_for_tuning_budget(hard_timeout)

  h0_budget <- min(r, max(0.05, 0.45 * hard_timeout))
  ev0 <- eval_theta_timebudget(theta = theta, time_budget_sec = h0_budget,
                               alpha = alpha, scenario = "H0", rH1 = rH1,
                               min_runs = min_h0_runs, deadline = deadline,
                               max_n_total = max_n_total, use_fast_ks = use_fast_ks)

  valid_type1 <- is.finite(ev0$rate) && ev0$n_run >= min_h0_runs && ev0$rate <= alpha
  remaining <- max(0, deadline - now_elapsed())

  ev1 <- if (valid_type1 && remaining > 0.05) {
    eval_theta_timebudget(theta = theta, time_budget_sec = min(r, remaining),
                          alpha = alpha, scenario = "H1", rH1 = rH1,
                          min_runs = min_h1_runs, deadline = deadline,
                          max_n_total = max_n_total, use_fast_ks = use_fast_ks)
  } else {
    list(rate = NA_real_, n_run = 0L, n_rej = 0L,
         avg_cpu_per_run = NA_real_, elapsed = 0,
         enough_runs = FALSE, timed_out = TRUE,
         n_aborted = 0L, abort_reasons = "")
  }

  enough_h1 <- valid_type1 && is.finite(ev1$rate) && ev1$n_run >= min_h1_runs
  power_hat <- if (enough_h1) ev1$rate else NA_real_
  loss_hat <- if (valid_type1 && enough_h1 && is.finite(power_hat)) 1 - power_hat else Inf
  timed <- isTRUE(ev0$timed_out) || (valid_type1 && isTRUE(ev1$timed_out)) ||
    !ev0$enough_runs || (valid_type1 && !enough_h1)

  list(valid = valid_type1 && ev0$enough_runs,
       type1_hat = ev0$rate,
       type1_n = ev0$n_run,
       cpu0_hat = ev0$avg_cpu_per_run,
       power_hat = power_hat,
       power_n = ev1$n_run,
       power_rej = ev1$n_rej,
       loss = loss_hat,
       timed_out = timed,
       h0_timed_out = isTRUE(ev0$timed_out),
       h1_timed_out = isTRUE(ev1$timed_out))
}

# -----------------------------------------
# Uniform allocation
# -----------------------------------------
tune_uniform_timebudgeted <- function(thetas,
                                      B_tune_sec,
                                      r_eval = 0.25,
                                      alpha = 0.05,
                                      rH1 = rH1_shift005,
                                      seed = 1,
                                      shuffle_initial = TRUE) {
  set.seed(seed)
  M <- length(thetas)
  order_idx <- seq_len(M)
  if (shuffle_initial) order_idx <- sample(order_idx, M, replace = FALSE)

  t0 <- now_elapsed()
  n_theta_eval <- 0L
  n_timed_out <- 0L
  seen_idx <- integer(0)
  cycle_id <- 0L
  best_idx <- NA_integer_
  best_loss <- Inf
  best_theta <- thetas[[order_idx[1]]]

  while (remaining_budget(t0, B_tune_sec) > 0.05) {
    cycle_id <- cycle_id + 1L
    for (i in order_idx) {
      rem <- remaining_budget(t0, B_tune_sec)
      if (rem <= 0.05) break
      ev <- evaluate_theta_valid(theta = thetas[[i]], r = min(r_eval, rem),
                                 alpha = alpha, rH1 = rH1,
                                 hard_timeout = min(2.5, rem))
      n_theta_eval <- n_theta_eval + 1L
      seen_idx <- c(seen_idx, i)
      if (isTRUE(ev$timed_out)) n_timed_out <- n_timed_out + 1L
      if (ev$valid && is.finite(ev$loss) && ev$loss < best_loss) {
        best_loss <- ev$loss
        best_idx <- i
        best_theta <- thetas[[i]]
      }
      if (remaining_budget(t0, B_tune_sec) <= 0.05) break
    }
  }

  if (!is.finite(best_loss)) {
    best_idx <- order_idx[1]
    best_theta <- thetas[[best_idx]]
  }

  list(best_theta = best_theta, best_idx = best_idx, best_loss = best_loss,
       method = "Uniform", tune_elapsed = as.numeric(now_elapsed() - t0),
       n_theta_eval = n_theta_eval, n_theta_eval_unique = length(unique(seen_idx)),
       n_timed_out = n_timed_out, r_eval = r_eval, cycles_completed = cycle_id)
}

# -----------------------------------------
# Random search
# -----------------------------------------
tune_random_search_valid <- function(thetas,
                                     B_tune_sec,
                                     r_eval = 0.25,
                                     alpha = 0.05,
                                     rH1 = rH1_shift005,
                                     seed = 1) {
  set.seed(seed)
  M <- length(thetas)
  t0 <- now_elapsed()
  best_idx <- NA_integer_
  best_loss <- Inf
  best_theta <- thetas[[1]]
  n_eval <- 0L
  n_timed_out <- 0L
  seen_idx <- integer(0)

  while (remaining_budget(t0, B_tune_sec) > 0.05) {
    rem <- remaining_budget(t0, B_tune_sec)
    i <- sample.int(M, 1L)
    seen_idx <- c(seen_idx, i)
    ev <- evaluate_theta_valid(theta = thetas[[i]], r = min(r_eval, rem),
                               alpha = alpha, rH1 = rH1,
                               hard_timeout = min(2.5, rem))
    n_eval <- n_eval + 1L
    if (isTRUE(ev$timed_out)) n_timed_out <- n_timed_out + 1L
    if (ev$valid && is.finite(ev$loss) && ev$loss < best_loss) {
      best_loss <- ev$loss
      best_idx <- i
      best_theta <- thetas[[i]]
    }
  }

  if (!is.finite(best_loss)) {
    best_idx <- 1L
    best_theta <- thetas[[best_idx]]
  }

  list(best_theta = best_theta, best_idx = best_idx, best_loss = best_loss,
       method = "Random", tune_elapsed = as.numeric(now_elapsed() - t0),
       n_theta_eval = n_eval, n_theta_eval_unique = length(unique(seen_idx)),
       n_timed_out = n_timed_out, r_eval = r_eval)
}

# -----------------------------------------
# Simulated annealing on fixed candidate-set graph
# -----------------------------------------
tune_sa_valid_graph <- function(thetas,
                                B_tune_sec,
                                init_idx = NULL,
                                r_eval = 0.25,
                                T0 = 0.20,
                                cooling = 0.95,
                                alpha = 0.05,
                                rH1 = rH1_shift005,
                                seed = 1,
                                k_neighbors = 8L,
                                neighbors = NULL,
                                method_label = "SA") {
  set.seed(seed)
  M <- length(thetas)
  stopifnot(M >= 1L)
  if (is.null(neighbors)) neighbors <- build_knn_neighbors(thetas, k_neighbors = k_neighbors)

  t0 <- now_elapsed()
  cur_idx <- if (!is.null(init_idx)) as.integer(init_idx) else sample.int(M, 1L)
  if (!is.finite(cur_idx) || cur_idx < 1L || cur_idx > M) cur_idx <- sample.int(M, 1L)
  init_idx_original <- cur_idx
  cur_theta <- thetas[[cur_idx]]

  rem <- remaining_budget(t0, B_tune_sec)
  cur_ev <- evaluate_theta_valid(theta = cur_theta, r = min(r_eval, rem),
                                 alpha = alpha, rH1 = rH1,
                                 hard_timeout = min(2.5, rem))
  cur_loss <- if (cur_ev$valid) cur_ev$loss else Inf
  best_theta <- cur_theta
  best_idx <- cur_idx
  best_loss <- cur_loss
  iter <- 0L
  n_timed_out <- as.integer(isTRUE(cur_ev$timed_out))
  T <- T0

  repeat {
    rem <- remaining_budget(t0, B_tune_sec)
    if (rem <= 0.05) break
    iter <- iter + 1L
    prop_idx <- propose_neighbor_graph(cur_idx, neighbors)
    prop_theta <- thetas[[prop_idx]]
    prop_ev <- evaluate_theta_valid(theta = prop_theta, r = min(r_eval, rem),
                                    alpha = alpha, rH1 = rH1,
                                    hard_timeout = min(2.5, rem))
    if (isTRUE(prop_ev$timed_out)) n_timed_out <- n_timed_out + 1L
    prop_loss <- if (prop_ev$valid) prop_ev$loss else Inf

    accept <- FALSE
    if (!is.finite(cur_loss) && is.finite(prop_loss)) accept <- TRUE
    else if (!is.finite(cur_loss) && !is.finite(prop_loss)) accept <- TRUE
    else if (is.finite(prop_loss) && is.finite(cur_loss)) {
      dE <- prop_loss - cur_loss
      accept <- (dE <= 0) || (runif(1) < exp(-dE / max(1e-8, T)))
    }

    if (accept) {
      cur_idx <- prop_idx
      cur_theta <- prop_theta
      cur_loss <- prop_loss
    }
    if (is.finite(prop_loss) && prop_loss < best_loss) {
      best_loss <- prop_loss
      best_theta <- prop_theta
      best_idx <- prop_idx
    }
    T <- T * cooling
  }

  list(best_theta = best_theta, best_idx = best_idx, best_loss = best_loss,
       method = method_label, tune_elapsed = as.numeric(now_elapsed() - t0),
       n_theta_eval = iter + 1L, n_timed_out = n_timed_out,
       r_eval = r_eval, T0 = T0, cooling = cooling,
       init_idx = init_idx_original, k_neighbors = k_neighbors)
}

tune_sa_valid <- function(thetas,
                          levels = NULL,
                          B_tune_sec,
                          r_eval = 0.25,
                          T0 = 0.20,
                          cooling = 0.95,
                          alpha = 0.05,
                          rH1 = rH1_shift005,
                          seed = 1) {
  tune_sa_valid_graph(thetas = thetas, B_tune_sec = B_tune_sec,
                      init_idx = NULL, r_eval = r_eval, T0 = T0,
                      cooling = cooling, alpha = alpha, rH1 = rH1,
                      seed = seed, method_label = "SA")
}

tune_sa_valid_from_init <- function(thetas,
                                    levels = NULL,
                                    init_idx = NULL,
                                    init_theta = NULL,
                                    B_tune_sec,
                                    r_eval = 0.25,
                                    T0 = 0.20,
                                    cooling = 0.95,
                                    alpha = 0.05,
                                    rH1 = rH1_shift005,
                                    seed = 1) {
  if (is.null(init_idx) && !is.null(init_theta)) {
    lookup <- build_theta_lookup(thetas)
    init_idx <- unname(lookup[theta_id_string(init_theta)])
    if (!is.finite(init_idx)) stop("init_theta is not found in the supplied candidate set.")
  }
  tune_sa_valid_graph(thetas = thetas, B_tune_sec = B_tune_sec,
                      init_idx = init_idx, r_eval = r_eval, T0 = T0,
                      cooling = cooling, alpha = alpha, rH1 = rH1,
                      seed = seed, method_label = "SA_from_init")
}

# -----------------------------------------
# Successive Halving
# -----------------------------------------
successive_halving_timebudgeted <- function(thetas,
                                            B_tune_sec,
                                            eta = 2,
                                            R = NULL,
                                            round_budgets = NULL,
                                            fidelity_ladder = NULL,
                                            alpha = 0.05,
                                            rH1 = rH1_shift005,
                                            lcb_conf = 0.90,
                                            seed = 1,
                                            shuffle_each_round = TRUE,
                                            require_complete_round = TRUE) {
  set.seed(seed)
  stopifnot(eta > 1)
  M <- length(thetas)
  stopifnot(M >= 2)

  if (is.null(R)) R <- max(1L, min(ceiling(log(M, base = eta)), 6L))
  if (is.null(round_budgets)) {
    round_budgets <- rep(B_tune_sec / R, R)
  } else {
    stopifnot(length(round_budgets) == R)
    if (sum(round_budgets) > B_tune_sec) stop("sum(round_budgets) must be <= B_tune_sec")
  }
  if (is.null(fidelity_ladder)) {
    r_max <- 0.60
    r1 <- r_max / eta^(R - 1)
    fidelity_ladder <- r1 * eta^(0:(R - 1))
  } else {
    stopifnot(length(fidelity_ladder) == R)
  }

  S <- seq_len(M)
  history <- list()
  best_idx <- NA_integer_
  best_loss <- Inf
  best_rank_loss <- Inf
  best_theta <- thetas[[1]]
  t0_outer <- now_elapsed()
  n_theta_eval_total <- 0L
  n_timed_out <- 0L

  for (r in seq_len(R)) {
    if (length(S) == 0L) break
    if (remaining_budget(t0_outer, B_tune_sec) <= 0.05) break
    S_round <- S
    if (shuffle_each_round) S_round <- sample(S_round, length(S_round), replace = FALSE)

    rk <- fidelity_ladder[r]
    t0_round <- now_elapsed()
    rows <- list()
    evaluated_idx <- integer(0)

    for (pos in seq_along(S_round)) {
      i <- S_round[pos]
      rem_outer <- remaining_budget(t0_outer, B_tune_sec)
      rem_round <- max(0, round_budgets[r] - (now_elapsed() - t0_round))
      rem <- min(rem_outer, rem_round)
      if (rem <= 0.05) break

      # Increase fidelity slightly when the round has spare budget, while
      # keeping the originally intended ladder as a lower bound.
      candidates_left <- length(S_round) - pos + 1L
      rk_dynamic <- min(rem, max(rk, rem_round / max(1L, candidates_left)))

      ev <- evaluate_theta_valid(theta = thetas[[i]], r = min(rk_dynamic, rem),
                                 alpha = alpha, rH1 = rH1,
                                 hard_timeout = min(2.5, rem))
      n_theta_eval_total <- n_theta_eval_total + 1L
      evaluated_idx <- c(evaluated_idx, i)
      if (isTRUE(ev$timed_out)) n_timed_out <- n_timed_out + 1L

      power_lcb <- if (ev$valid && ev$power_n > 0) wilson_lower(ev$power_rej, ev$power_n, conf = lcb_conf) else 0
      rank_loss_i <- if (ev$valid) 1 - power_lcb else Inf
      score_loss_i <- if (ev$valid && is.finite(ev$power_hat)) ev$loss else Inf

      rows[[length(rows) + 1L]] <- data.frame(
        idx = i, valid = ev$valid, timed_out = ev$timed_out,
        type1_hat = ev$type1_hat, cpu0_hat = ev$cpu0_hat,
        power_hat = ev$power_hat, power_n = ev$power_n,
        power_rej = ev$power_rej, power_lcb = power_lcb,
        rank_loss = rank_loss_i, loss = score_loss_i,
        r_eval = rk_dynamic, stringsAsFactors = FALSE)

      if (is.finite(score_loss_i) && score_loss_i < best_loss) {
        best_loss <- score_loss_i
        best_rank_loss <- rank_loss_i
        best_idx <- i
        best_theta <- thetas[[i]]
      }
    }

    tab <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
    truncated <- length(evaluated_idx) < length(S)

    if (nrow(tab) == 0L) {
      history[[r]] <- list(round = r, r_eval = rk, round_budget = round_budgets[r],
                           table = tab, survivors = S, n_eval = 0L, truncated = TRUE)
      break
    }

    tab <- tab[order(tab$rank_loss, -tab$power_n), , drop = FALSE]

    if (truncated && require_complete_round) {
      history[[r]] <- list(round = r, r_eval = rk, round_budget = round_budgets[r],
                           table = tab, survivors = S, n_eval = nrow(tab), truncated = TRUE)
      break
    }

    keep <- max(1L, ceiling(nrow(tab) / eta))
    S_next <- tab$idx[seq_len(keep)]
    history[[r]] <- list(round = r, r_eval = rk, round_budget = round_budgets[r],
                         table = tab, survivors = S_next, n_eval = nrow(tab), truncated = FALSE)
    S <- S_next
    if (length(S) <= 1L) break
  }

  if (is.na(best_idx) && length(S) >= 1L) {
    best_idx <- S[1]
    best_theta <- thetas[[best_idx]]
  }

  list(best_theta = best_theta, best_idx = best_idx,
       best_loss = best_loss, best_rank_loss = best_rank_loss,
       history = history, K_rounds = length(history), method = "SH",
       tune_elapsed = as.numeric(now_elapsed() - t0_outer),
       n_theta_eval = n_theta_eval_total, n_timed_out = n_timed_out,
       eta = eta, R = R, round_budgets = round_budgets,
       fidelity_ladder = fidelity_ladder, lcb_conf = lcb_conf)
}

# -----------------------------------------
# Hybrid: SH first, then SA refinement
# -----------------------------------------
tune_hybrid_sh_sa <- function(thetas,
                              levels = NULL,
                              B_tune_sec,
                              lambda_sh = 0.50,
                              sh_eta = 2,
                              sh_lcb_conf = 0.90,
                              r_eval = 0.25,
                              sa_T0 = 0.20,
                              sa_cooling = 0.95,
                              alpha = 0.05,
                              rH1 = rH1_shift005,
                              seed = 1,
                              restrict_to_sh_survivors = FALSE,
                              k_neighbors = 8L) {
  stopifnot(lambda_sh > 0 && lambda_sh < 1)
  B_sh <- lambda_sh * B_tune_sec
  t0_total <- now_elapsed()

  out_sh <- successive_halving_timebudgeted(
    thetas = thetas, B_tune_sec = B_sh, eta = sh_eta,
    alpha = alpha, rH1 = rH1, lcb_conf = sh_lcb_conf, seed = seed)

  init_obj <- choose_hybrid_init_from_sh(out_sh, thetas)
  B_sa <- remaining_budget(t0_total, B_tune_sec)

  if (B_sa <= 0.05) {
    out_sa <- list(best_theta = init_obj$init_theta, best_idx = init_obj$init_idx,
                   best_loss = Inf, tune_elapsed = 0, n_theta_eval = 0L,
                   n_timed_out = 0L, init_idx = init_obj$init_idx)
  } else if (restrict_to_sh_survivors) {
    survivors <- last_survivors_from_sh(out_sh)
    survivors <- unique(survivors[is.finite(survivors) & survivors >= 1L & survivors <= length(thetas)])
    if (length(survivors) < 2L) survivors <- unique(c(init_obj$init_idx, seq_len(min(2L, length(thetas)))))
    thetas_sa <- thetas[survivors]
    init_idx_sa <- match(init_obj$init_idx, survivors)
    if (!is.finite(init_idx_sa)) init_idx_sa <- 1L

    out_sa_local <- tune_sa_valid_graph(
      thetas = thetas_sa, B_tune_sec = B_sa, init_idx = init_idx_sa,
      r_eval = r_eval, T0 = sa_T0, cooling = sa_cooling,
      alpha = alpha, rH1 = rH1, seed = seed,
      k_neighbors = min(k_neighbors, max(1L, length(thetas_sa) - 1L)),
      method_label = "SA_from_SH_restricted")

    out_sa <- out_sa_local
    out_sa$best_idx_local <- out_sa_local$best_idx
    out_sa$best_idx <- survivors[out_sa_local$best_idx]
    out_sa$best_theta <- thetas[[out_sa$best_idx]]
    out_sa$init_idx_local <- init_idx_sa
    out_sa$init_idx <- init_obj$init_idx
    out_sa$survivor_set_size <- length(survivors)
  } else {
    out_sa <- tune_sa_valid_graph(
      thetas = thetas, B_tune_sec = B_sa, init_idx = init_obj$init_idx,
      r_eval = r_eval, T0 = sa_T0, cooling = sa_cooling,
      alpha = alpha, rH1 = rH1, seed = seed,
      k_neighbors = k_neighbors, method_label = "SA_from_SH")
  }

  best_idx <- out_sh$best_idx
  best_theta <- out_sh$best_theta
  best_loss <- out_sh$best_loss
  if (is.finite(out_sa$best_loss) && (!is.finite(best_loss) || out_sa$best_loss < best_loss)) {
    best_idx <- out_sa$best_idx
    best_theta <- out_sa$best_theta
    best_loss <- out_sa$best_loss
  }

  list(best_theta = best_theta, best_idx = best_idx, best_loss = best_loss,
       method = if (restrict_to_sh_survivors) "Hybrid_restricted" else "Hybrid",
       tune_elapsed = as.numeric(now_elapsed() - t0_total),
       n_theta_eval = out_sh$n_theta_eval + out_sa$n_theta_eval,
       n_timed_out = out_sh$n_timed_out + out_sa$n_timed_out,
       lambda_sh = lambda_sh, B_sh = B_sh, B_sa = B_sa,
       sh_stage = out_sh, sa_stage = out_sa,
       sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
       sa_T0 = sa_T0, sa_cooling = sa_cooling,
       r_eval = r_eval, sa_init_idx = init_obj$init_idx,
       restrict_to_sh_survivors = restrict_to_sh_survivors)
}

# -----------------------------------------
# Final reporting row
# -----------------------------------------
final_row_from_tuner <- function(out, method_name, seed, theta_set_id, theta_set_size,
                                 B_tune_sec, alpha_in = 0.05,
                                 N0_final = 500, N1_final = 500) {
  set.seed(stable_seed(theta_set_id, B_tune_sec, method_name, seed, "H0_final"))
  fe0 <- final_eval_fixedN(out$best_theta, N = N0_final, alpha = alpha_in, scenario = "H0")
  set.seed(stable_seed(theta_set_id, B_tune_sec, method_name, seed, "H1_final"))
  fe1 <- final_eval_fixedN(out$best_theta, N = N1_final, alpha = alpha_in, scenario = "H1")

  valid_flag <- as.integer(is.finite(fe0$rate) && fe0$rate <= alpha_in)
  tune_elapsed <- as.numeric(out$tune_elapsed)

  data.frame(
    seed = seed, theta_set_id = theta_set_id,
    theta_set_size = theta_set_size, B_tune_sec = B_tune_sec,
    method = method_name,
    best_theta_n0 = out$best_theta$n0,
    best_theta_gamma = out$best_theta$gamma,
    best_theta_K = out$best_theta$K,
    best_theta_rho = out$best_theta$rho,
    best_theta_u = out$best_theta$u,
    tune_elapsed_sec = tune_elapsed,
    budget_overrun_sec = pmax(tune_elapsed - B_tune_sec, 0),
    budget_overrun_ratio = ifelse(B_tune_sec > 0, tune_elapsed / B_tune_sec, NA_real_),
    n_theta_eval = out$n_theta_eval,
    n_theta_eval_unique = if (!is.null(out$n_theta_eval_unique)) out$n_theta_eval_unique else NA_integer_,
    n_timed_out = out$n_timed_out,
    lambda_sh = if (!is.null(out$lambda_sh)) out$lambda_sh else NA_real_,
    sh_eta = if (!is.null(out$sh_eta)) out$sh_eta else NA_real_,
    sh_lcb_conf = if (!is.null(out$sh_lcb_conf)) out$sh_lcb_conf else NA_real_,
    sa_T0 = if (!is.null(out$sa_T0)) out$sa_T0 else NA_real_,
    sa_cooling = if (!is.null(out$sa_cooling)) out$sa_cooling else NA_real_,
    r_eval = if (!is.null(out$r_eval)) out$r_eval else NA_real_,
    final_N0 = N0_final,
    final_type1_hat = fe0$rate,
    final_type1_ci_low = fe0$ci_low,
    final_type1_ci_high = fe0$ci_high,
    final_N1 = N1_final,
    final_power_hat = fe1$rate,
    final_power_ci_low = fe1$ci_low,
    final_power_ci_high = fe1$ci_high,
    final_c0_mean = fe0$cpu_mean,
    final_c0_median = fe0$cpu_median,
    final_c0_q90 = fe0$cpu_q90,
    final_c0_q95 = fe0$cpu_q95,
    valid_flag = valid_flag,
    valid_power = fe1$rate * valid_flag,
    type1_excess = pmax(fe0$rate - alpha_in, 0),
    stringsAsFactors = FALSE)
}

# -----------------------------------------
# Main method comparison for one regime / one seed
# -----------------------------------------
run_one_comparison_with_hybrid <- function(seed = 1,
                                           thetas_in,
                                           levels_in = make_levels_from_thetas(thetas_in),
                                           theta_set_id = NA_character_,
                                           theta_set_size = length(thetas_in),
                                           B_tune_in = 160,
                                           alpha_in = 0.05,
                                           N0_final = 500,
                                           N1_final = 500,
                                           r_eval_baseline = 0.25,
                                           sh_eta = 2,
                                           sh_lcb_conf = 0.90,
                                           sa_T0 = 0.20,
                                           sa_cooling = 0.95,
                                           lambda_sh = 0.50) {
  set.seed(seed)

  out_sh <- successive_halving_timebudgeted(
    thetas = thetas_in, B_tune_sec = B_tune_in, eta = sh_eta,
    alpha = alpha_in, rH1 = rH1_shift005, lcb_conf = sh_lcb_conf, seed = seed)

  out_uni <- tune_uniform_timebudgeted(
    thetas = thetas_in, B_tune_sec = B_tune_in, r_eval = r_eval_baseline,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed, shuffle_initial = TRUE)

  out_rs <- tune_random_search_valid(
    thetas = thetas_in, B_tune_sec = B_tune_in, r_eval = r_eval_baseline,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed)

  out_sa <- tune_sa_valid(
    thetas = thetas_in, levels = levels_in, B_tune_sec = B_tune_in,
    r_eval = r_eval_baseline, T0 = sa_T0, cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed)

  out_hybrid <- tune_hybrid_sh_sa(
    thetas = thetas_in, levels = levels_in, B_tune_sec = B_tune_in,
    lambda_sh = lambda_sh, sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
    r_eval = r_eval_baseline, sa_T0 = sa_T0, sa_cooling = sa_cooling,
    alpha = alpha_in, rH1 = rH1_shift005, seed = seed)

  outs <- list(SH = out_sh, Uniform = out_uni, Random = out_rs,
               SA = out_sa, Hybrid = out_hybrid)

  bind_rows(lapply(names(outs), function(meth) {
    final_row_from_tuner(outs[[meth]], method_name = meth,
                         seed = seed, theta_set_id = theta_set_id,
                         theta_set_size = theta_set_size,
                         B_tune_sec = B_tune_in,
                         alpha_in = alpha_in,
                         N0_final = N0_final, N1_final = N1_final)
  }))
}

run_experiment_grid_with_hybrid <- function(
    all_thetas = full_thetas,
    candidate_sizes = c(64, 128, 256, 512),
    budget_grid = c(40, 80, 160, 320),
    seeds = 1:5,
    alpha_in = 0.05,
    N0_final = 500,
    N1_final = 500,
    r_eval_baseline = 0.25,
    sh_eta = 2,
    sh_lcb_conf = 0.90,
    sa_T0 = 0.20,
    sa_cooling = 0.95,
    lambda_sh = 0.50,
    output_csv = "results_hybrid_main.csv") {

  stopifnot(max(candidate_sizes) <= length(all_thetas))
  all_rows <- list()
  row_id <- 0L

  for (M in candidate_sizes) {
    for (s in seeds) {
      subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
      theta_set_id <- paste0("M", M, "_seed", s)
      for (B_tune in budget_grid) {
        row_id <- row_id + 1L
        cat("Running regime: |Theta| =", M,
            ", B_tune =", B_tune, ", seed =", s, "\n")
        all_rows[[row_id]] <- run_one_comparison_with_hybrid(
          seed = s, thetas_in = subset_obj$thetas, levels_in = subset_obj$levels,
          theta_set_id = theta_set_id, theta_set_size = M, B_tune_in = B_tune,
          alpha_in = alpha_in, N0_final = N0_final, N1_final = N1_final,
          r_eval_baseline = r_eval_baseline, sh_eta = sh_eta,
          sh_lcb_conf = sh_lcb_conf, sa_T0 = sa_T0,
          sa_cooling = sa_cooling, lambda_sh = lambda_sh)
      }
    }
  }

  out_df <- bind_rows(all_rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

# -----------------------------------------
# Sensitivity helpers retained for earlier chapters / appendices
# -----------------------------------------
run_sh_sensitivity <- function(
    all_thetas = full_thetas,
    candidate_sizes = c(64, 512),
    budget_grid = c(40, 320),
    seeds = 1:5,
    eta_grid = c(2, 3, 4),
    alpha_in = 0.05,
    N0_final = 300,
    N1_final = 300,
    sh_lcb_conf = 0.90,
    output_csv = "results_sh_sensitivity.csv") {

  rows <- list(); id <- 0L
  for (M in candidate_sizes) for (s in seeds) {
    subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
    theta_set_id <- paste0("SHsens_M", M, "_seed", s)
    for (B_tune in budget_grid) for (eta_val in eta_grid) {
      cat("SH sensitivity: |Theta| =", M, ", B =", B_tune,
          ", seed =", s, ", eta =", eta_val, "\n")
      out <- successive_halving_timebudgeted(
        thetas = subset_obj$thetas, B_tune_sec = B_tune,
        eta = eta_val, alpha = alpha_in, rH1 = rH1_shift005,
        lcb_conf = sh_lcb_conf, seed = s)
      id <- id + 1L
      rows[[id]] <- final_row_from_tuner(out, "SH", s, theta_set_id, M,
                                         B_tune, alpha_in, N0_final, N1_final)
      rows[[id]]$sh_eta <- eta_val
      rows[[id]]$sh_lcb_conf <- sh_lcb_conf
    }
  }
  out_df <- bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

run_sa_sensitivity <- function(
    all_thetas = full_thetas,
    candidate_sizes = c(64, 512),
    budget_grid = c(40, 320),
    seeds = 1:5,
    T0_grid = c(0.10, 0.20, 0.50),
    cooling_grid = c(0.90, 0.95),
    alpha_in = 0.05,
    N0_final = 300,
    N1_final = 300,
    r_eval_baseline = 0.25,
    output_csv = "results_sa_sensitivity.csv") {

  rows <- list(); id <- 0L
  for (M in candidate_sizes) for (s in seeds) {
    subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
    theta_set_id <- paste0("SAsens_M", M, "_seed", s)
    for (B_tune in budget_grid) for (T0_val in T0_grid) for (cool_val in cooling_grid) {
      cat("SA sensitivity: |Theta| =", M, ", B =", B_tune,
          ", seed =", s, ", T0 =", T0_val, ", cooling =", cool_val, "\n")
      out <- tune_sa_valid(
        thetas = subset_obj$thetas, levels = subset_obj$levels,
        B_tune_sec = B_tune, r_eval = r_eval_baseline,
        T0 = T0_val, cooling = cool_val,
        alpha = alpha_in, rH1 = rH1_shift005, seed = s)
      id <- id + 1L
      rows[[id]] <- final_row_from_tuner(out, "SA", s, theta_set_id, M,
                                         B_tune, alpha_in, N0_final, N1_final)
      rows[[id]]$sa_T0 <- T0_val
      rows[[id]]$sa_cooling <- cool_val
    }
  }
  out_df <- bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

# -----------------------------------------
# Hybrid internal tuning retained for compatibility
# -----------------------------------------
make_hybrid_param_grid <- function(
    lambda_sh_grid   = c(0.30, 0.50, 0.70),
    sh_eta_grid      = c(2, 3),
    sa_T0_grid       = c(0.10, 0.20, 0.50),
    sa_cooling_grid  = c(0.90, 0.95),
    sh_lcb_conf_grid = c(0.90),
    r_eval_grid      = c(0.25)) {
  expand.grid(lambda_sh = lambda_sh_grid, sh_eta = sh_eta_grid,
              sa_T0 = sa_T0_grid, sa_cooling = sa_cooling_grid,
              sh_lcb_conf = sh_lcb_conf_grid, r_eval = r_eval_grid,
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

run_one_hybrid_setting <- function(seed = 1,
                                   thetas_in,
                                   levels_in = make_levels_from_thetas(thetas_in),
                                   theta_set_id = NA_character_,
                                   theta_set_size = length(thetas_in),
                                   B_tune_in = 160,
                                   alpha_in = 0.05,
                                   N0_final = 300,
                                   N1_final = 300,
                                   lambda_sh = 0.50,
                                   sh_eta = 2,
                                   sh_lcb_conf = 0.90,
                                   r_eval = 0.25,
                                   sa_T0 = 0.20,
                                   sa_cooling = 0.95) {
  out <- tune_hybrid_sh_sa(thetas = thetas_in, levels = levels_in,
                           B_tune_sec = B_tune_in, lambda_sh = lambda_sh,
                           sh_eta = sh_eta, sh_lcb_conf = sh_lcb_conf,
                           r_eval = r_eval, sa_T0 = sa_T0,
                           sa_cooling = sa_cooling, alpha = alpha_in,
                           rH1 = rH1_shift005, seed = seed)
  row <- final_row_from_tuner(out, "Hybrid", seed, theta_set_id,
                              theta_set_size, B_tune_in,
                              alpha_in, N0_final, N1_final)
  row$sh_stage_n_eval <- out$sh_stage$n_theta_eval
  row$sa_stage_n_eval <- out$sa_stage$n_theta_eval
  row$sa_init_idx <- out$sa_init_idx
  row
}

run_hybrid_internal_tuning <- function(
    all_thetas = full_thetas,
    candidate_sizes = c(128, 512),
    budget_grid = c(80, 160, 320),
    seeds = 1:3,
    param_grid = make_hybrid_param_grid(),
    alpha_in = 0.05,
    N0_final = 300,
    N1_final = 300,
    output_csv = "results_hybrid_internal_tuning.csv") {
  rows <- list(); id <- 0L
  cat("Total hybrid tuning jobs =",
      nrow(param_grid) * length(candidate_sizes) * length(budget_grid) * length(seeds), "\n")
  for (g in seq_len(nrow(param_grid))) {
    pars <- param_grid[g, ]
    for (M in candidate_sizes) for (s in seeds) {
      subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
      theta_set_id <- paste0("TUNE_M", M, "_seed", s)
      for (B_tune in budget_grid) {
        cat("Hybrid tuning: setting =", g, ", |Theta| =", M,
            ", B =", B_tune, ", seed =", s, "\n")
        id <- id + 1L
        rows[[id]] <- run_one_hybrid_setting(
          seed = s, thetas_in = subset_obj$thetas, levels_in = subset_obj$levels,
          theta_set_id = theta_set_id, theta_set_size = M, B_tune_in = B_tune,
          alpha_in = alpha_in, N0_final = N0_final, N1_final = N1_final,
          lambda_sh = pars$lambda_sh, sh_eta = pars$sh_eta,
          sh_lcb_conf = pars$sh_lcb_conf, r_eval = pars$r_eval,
          sa_T0 = pars$sa_T0, sa_cooling = pars$sa_cooling)
      }
    }
  }
  out_df <- bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

summarise_hybrid_internal_tuning <- function(df) {
  df %>%
    group_by(lambda_sh, sh_eta, sh_lcb_conf, r_eval, sa_T0, sa_cooling) %>%
    summarise(n_runs = n(),
              mean_power = mean(final_power_hat, na.rm = TRUE),
              sd_power = sd(final_power_hat, na.rm = TRUE),
              mean_valid_rate = mean(valid_flag, na.rm = TRUE),
              mean_valid_power = mean(valid_power, na.rm = TRUE),
              mean_type1 = mean(final_type1_hat, na.rm = TRUE),
              mean_type1_excess = mean(type1_excess, na.rm = TRUE),
              mean_c0 = mean(final_c0_mean, na.rm = TRUE),
              mean_n_eval = mean(n_theta_eval, na.rm = TRUE),
              mean_n_timed_out = mean(n_timed_out, na.rm = TRUE),
              mean_sh_stage_n_eval = mean(sh_stage_n_eval, na.rm = TRUE),
              mean_sa_stage_n_eval = mean(sa_stage_n_eval, na.rm = TRUE),
              n_valid = sum(valid_flag, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_valid_power), desc(mean_valid_rate),
            mean_type1_excess, desc(mean_power), mean_c0)
}

select_top_hybrid_settings <- function(summary_df, top_k = 3) {
  summary_df[seq_len(min(top_k, nrow(summary_df))), , drop = FALSE]
}

validate_one_hybrid_setting_vs_baselines <- function(
    setting_row,
    all_thetas = full_thetas,
    candidate_sizes = c(64, 128, 256, 512),
    budget_grid = c(40, 80, 160, 320),
    seeds = 1:5,
    alpha_in = 0.05,
    N0_final = 500,
    N1_final = 500,
    output_csv = "results_best_hybrid_vs_baselines.csv") {
  stopifnot(nrow(setting_row) == 1)
  out_df <- run_experiment_grid_with_hybrid(
    all_thetas = all_thetas, candidate_sizes = candidate_sizes,
    budget_grid = budget_grid, seeds = seeds, alpha_in = alpha_in,
    N0_final = N0_final, N1_final = N1_final,
    r_eval_baseline = setting_row$r_eval, sh_eta = setting_row$sh_eta,
    sh_lcb_conf = setting_row$sh_lcb_conf, sa_T0 = setting_row$sa_T0,
    sa_cooling = setting_row$sa_cooling, lambda_sh = setting_row$lambda_sh,
    output_csv = output_csv)
  out_df
}

validate_top_hybrid_settings_only <- function(
    selected_settings,
    all_thetas = full_thetas,
    candidate_sizes = c(64, 128, 256, 512),
    budget_grid = c(40, 80, 160, 320),
    seeds = 1:5,
    alpha_in = 0.05,
    N0_final = 500,
    N1_final = 500,
    output_csv = "results_top_hybrid_settings_only.csv") {
  rows <- list(); id <- 0L
  for (k in seq_len(nrow(selected_settings))) {
    setting_row <- selected_settings[k, , drop = FALSE]
    for (M in candidate_sizes) for (s in seeds) {
      subset_obj <- sample_candidate_subset(all_thetas, M = M, seed = s)
      theta_set_id <- paste0("VALTOP", k, "_M", M, "_seed", s)
      for (B_tune in budget_grid) {
        id <- id + 1L
        rows[[id]] <- run_one_hybrid_setting(
          seed = s, thetas_in = subset_obj$thetas, levels_in = subset_obj$levels,
          theta_set_id = theta_set_id, theta_set_size = M, B_tune_in = B_tune,
          alpha_in = alpha_in, N0_final = N0_final, N1_final = N1_final,
          lambda_sh = setting_row$lambda_sh, sh_eta = setting_row$sh_eta,
          sh_lcb_conf = setting_row$sh_lcb_conf, r_eval = setting_row$r_eval,
          sa_T0 = setting_row$sa_T0, sa_cooling = setting_row$sa_cooling)
      }
    }
  }
  out_df <- bind_rows(rows)
  write.csv(out_df, output_csv, row.names = FALSE)
  out_df
}

# -----------------------------------------
# Budget diagnostics
# -----------------------------------------
budget_diagnostics <- function(df, tolerance_sec = 5) {
  if (!"budget_overrun_sec" %in% names(df)) {
    df$budget_overrun_sec <- pmax(df$tune_elapsed_sec - df$B_tune_sec, 0)
    df$budget_overrun_ratio <- df$tune_elapsed_sec / df$B_tune_sec
  }
  df %>%
    group_by(method) %>%
    summarise(n = n(),
              mean_elapsed = mean(tune_elapsed_sec, na.rm = TRUE),
              max_elapsed = max(tune_elapsed_sec, na.rm = TRUE),
              mean_spend_ratio = mean(budget_overrun_ratio, na.rm = TRUE),
              max_overrun = max(budget_overrun_sec, na.rm = TRUE),
              n_over_tolerance = sum(budget_overrun_sec > tolerance_sec, na.rm = TRUE),
              mean_n_timed_out = mean(n_timed_out, na.rm = TRUE),
              .groups = "drop")
}
