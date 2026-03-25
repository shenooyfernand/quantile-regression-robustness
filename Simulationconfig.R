# =============================================================================
# Simulation Configuration with Contamination & Heteroskedasticity
# =============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(future)
library(furrr)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)

# Reproducibility
set.seed(123)

# -----------------------------------------------------------------------------
# Simulation parameters
# -----------------------------------------------------------------------------
ns          = c(200, 500, 1000)

betas       = list(
  sparse       = c(5, 0, 0, 0, 0, 0, 0, 0),
  intermediate = c(3, 1.5, 0, 0, 2, 0, 0, 0),
  dense        = c(0, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0)
)

err_dists   = c("normal", "exp")
taus        = c(0.1, 0.25, 0.5, 0.75, 0.9)
losses      = c("pinball", "qhuber", "logcosh")

# Hyperparameter grids
qhuber_k_values  = c(0.05, 0.1, 0.2, 0.5, 1, 2, 3, 10)
logcosh_c_values = c(0.05, 0.1, 0.3, 0.5, 0.7, 1, 2, 10)

# Contamination parameters
contamination_levels = c(0, 0.05, 0.10, 0.20)
contam_scales        = c(10)

# Heteroskedasticity types
heterosked_types = c("none", "linear", "quadratic")

# =============================================================================
# MONTE CARLO SAMPLE SIZE JUSTIFICATION
# Following Morris et al. (2019) Section 5.3
# =============================================================================

# Define precision targets
mcse_target_coverage <- 0.01   # 1% for coverage (we don't have coverage yet)
mcse_target_bias     <- 0.01   # 0.01 for bias detection
mcse_target_conv     <- 0.01   # 1% for convergence rate

# Expected variability - conservative estimates
# Based on preliminary theory for n=200, tau=0.5
sd_beta_hat_est <- 0.25  # Conservative upper bound for SD(β̂)

# Calculate required n_rep for each criterion
cat("\n", rep("=", 70), "\n", sep = "")
cat("  MONTE CARLO SAMPLE SIZE JUSTIFICATION\n")
cat("  Following Morris et al. (2019) Section 5.3\n")
cat(rep("=", 70), "\n\n", sep = "")

# For bias precision
n_rep_bias <- ceiling((sd_beta_hat_est / mcse_target_bias)^2)
cat("For BIAS precision:\n")
cat("  Estimated SD(β̂):        ", sd_beta_hat_est, "\n")
cat("  Target MC SE:            ", mcse_target_bias, "\n")
cat("  Formula: n_rep = (SD(β̂) / target_MCSE)²\n")
cat("  Required n_rep:          ", n_rep_bias, "\n\n")

# For convergence rate precision (worst case: p = 0.5)
n_rep_conv <- ceiling(0.25 / mcse_target_conv^2)
cat("For CONVERGENCE RATE precision:\n")
cat("  Target MC SE:            ", mcse_target_conv, 
    sprintf("(%.1f%%)", mcse_target_conv*100), "\n")
cat("  Worst case: p = 0.5\n")
cat("  Formula: n_rep = 0.25 / target_MCSE²\n")
cat("  Required n_rep:          ", n_rep_conv, "\n\n")

# Choose the maximum (most conservative)
n_rep_calculated <- max(n_rep_bias, n_rep_conv)
cat("Calculated minimum n_rep:", n_rep_calculated, "\n")

# Monte Carlo replications - JUSTIFIED!
n_rep <- 1000  # Using 1000 for computational feasibility

cat("\n")
if (n_rep >= n_rep_calculated) {
  cat("✓ SELECTED n_rep =", n_rep, "(exceeds minimum)\n")
} else {
  cat("⚠ SELECTED n_rep =", n_rep, "(below calculated minimum of", 
      n_rep_calculated, ")\n")
  cat("  This provides MCSE(bias) ≈", 
      sprintf("%.4f", sd_beta_hat_est / sqrt(n_rep)), "\n")
  cat("  This provides MCSE(conv) ≈", 
      sprintf("%.4f (%.2f%%)", sqrt(0.25/n_rep), sqrt(0.25/n_rep)*100), "\n")
}

cat("\nActual Monte Carlo SEs achieved with n_rep =", n_rep, ":\n")
achieved_mcse_bias <- sd_beta_hat_est / sqrt(n_rep)
achieved_mcse_conv <- sqrt(0.25 / n_rep)
cat("  MCSE(bias):              ", sprintf("%.5f", achieved_mcse_bias), "\n")
cat("  MCSE(conv) worst case:   ", sprintf("%.5f (%.2f%%)", 
                                           achieved_mcse_conv, 
                                           achieved_mcse_conv*100), "\n")

cat("\n", rep("=", 70), "\n\n", sep = "")

# Metrics to track
metrics = c(
  "bias", "mse", "mae",
  "computation_time", "convergence", "iterations"
)
# Note: coverage and ci_width removed (no SEs computed)

# -----------------------------------------------------------------------------
# Build scenario grid
# -----------------------------------------------------------------------------

scenario_grid = expand.grid(
  n             = ns,
  beta_name     = names(betas),
  err_dist      = err_dists,
  contamination = contamination_levels,
  heterosked    = heterosked_types,
  stringsAsFactors = FALSE
) |>
  mutate(
    scenario_id  = dplyr::row_number(),
    beta_true    = purrr::map(beta_name, ~ betas[[.x]]),
    n_nonzero    = purrr::map_int(beta_true, ~ sum(.x != 0)),
    contam_scale = 10
  )

# -----------------------------------------------------------------------------
# Create full experimental grid
# -----------------------------------------------------------------------------

# 1) Pinball: no hyperparameters
grid_pinball <- scenario_grid |>
  tidyr::crossing(
    tau  = taus,
    loss = "pinball"
  ) |>
  mutate(
    k_param = NA_real_,
    c_param = NA_real_
  )

# 2) QHuber: expand only k_param
grid_qhuber <- scenario_grid |>
  tidyr::crossing(
    tau     = taus,
    loss    = "qhuber",
    k_param = qhuber_k_values
  ) |>
  mutate(
    c_param = NA_real_
  )

# 3) Log-cosh: expand only c_param
grid_logcosh <- scenario_grid |>
  tidyr::crossing(
    tau     = taus,
    loss    = "logcosh",
    c_param = logcosh_c_values
  ) |>
  mutate(
    k_param = NA_real_
  )

# Combine everything
experiment_grid <- dplyr::bind_rows(
  grid_pinball,
  grid_qhuber,
  grid_logcosh
) |>
  mutate(
    experiment_id = dplyr::row_number()
  )

# -----------------------------------------------------------------------------
# Configuration metadata
# -----------------------------------------------------------------------------

config_info = list(
  timestamp         = Sys.time(),
  r_version         = R.version.string,
  seed              = 123,
  n_scenarios       = nrow(scenario_grid),
  n_experiments     = nrow(experiment_grid),
  total_fits        = nrow(experiment_grid) * n_rep,
  estimated_hours   = (nrow(experiment_grid) * n_rep * 0.5) / 3600,
  packages          = names(sessionInfo()$otherPkgs),
  # NEW: Monte Carlo precision info
  mc_precision = list(
    n_rep = n_rep,
    n_rep_calculated = n_rep_calculated,
    target_mcse_bias = mcse_target_bias,
    target_mcse_conv = mcse_target_conv,
    achieved_mcse_bias = achieved_mcse_bias,
    achieved_mcse_conv = achieved_mcse_conv,
    justification = "See Morris et al. (2019) Section 5.3"
  ),
  contamination_info = list(
    levels = contamination_levels,
    scales = contam_scales
  ),
  heterosked_info = list(
    types = heterosked_types
  )
)

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

cat("=== Enhanced Simulation Configuration ===\n")
cat("Scenarios:                ", nrow(scenario_grid), "\n")
cat("  - Sample sizes:         ", paste(ns, collapse = ", "), "\n")
cat("  - Beta patterns:        ", paste(names(betas), collapse = ", "), "\n")
cat("  - Error distributions:  ", paste(err_dists, collapse = ", "), "\n")
cat("  - Contamination levels: ", paste(contamination_levels, collapse = ", "), "\n")
cat("  - Heterosked types:     ", paste(heterosked_types, collapse = ", "), "\n")
cat("Tau levels:               ", paste(taus, collapse = ", "), "\n")
cat("Losses:                   ", paste(losses, collapse = ", "), "\n")
cat("k values for QHuber:      ", paste(qhuber_k_values, collapse = ", "), "\n")
cat("c values for Log-cosh:    ", paste(logcosh_c_values, collapse = ", "), "\n")
cat("Total experiments:        ", nrow(experiment_grid), "\n")
cat("Replications/experiment:  ", n_rep, "\n")
cat("Total fits:               ", format(config_info$total_fits, big.mark = ","), "\n")
cat("Estimated time (rough):   ", round(config_info$estimated_hours, 1), "hours\n")
cat("Parallel workers:         ", future::nbrOfWorkers(), "\n")
cat("=========================================\n\n")

# -----------------------------------------------------------------------------
# Save configuration
# -----------------------------------------------------------------------------

save(
  scenario_grid,
  experiment_grid,
  config_info,
  betas,
  ns,
  err_dists,
  taus,
  losses,
  qhuber_k_values,
  logcosh_c_values,
  contamination_levels,
  contam_scales,
  heterosked_types,
  n_rep,
  metrics,
  file = "simulation_config.RData"
)

cat("Configuration saved to 'simulation_config.RData'\n")