# =============================================================================
# Small Test Simulation Configuration - Pilot with New Features
# =============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(furrr)
library(future)

plan(multisession)
set.seed(123)

# -----------------------------------------------------------------------------
# Minimal parameter set for testing
# -----------------------------------------------------------------------------

ns        = c(200)
betas     = list(sparse = c(5, 0, 0, 0, 0, 0, 0, 0))
err_dists = c("normal", "exp")  # Test both error distributions
taus      = c(0.5, 0.9)         # Test median and upper quantile

# Loss functions
losses_pinball = "pinball"
losses_qhuber  = "qhuber"
losses_logcosh = "logcosh"

# Hyperparameters (reduced grid for testing)
qhuber_k_values   = c(0.5, 1)
logcosh_c_values  = c(0.5, 1)

# NEW: Contamination and heteroskedasticity
contamination_levels = c(0, 0.10)      # Clean and 10% contamination
contam_scales        = c(10)           # Fixed outlier scale
heterosked_types     = c("none", "linear")  # Homoskedastic and heteroskedastic

# Monte Carlo replications (small for quick test)
n_rep = 10

# -----------------------------------------------------------------------------
# Build base scenario grid
# -----------------------------------------------------------------------------

scenario_grid <- expand.grid(
  n             = ns,
  beta_name     = names(betas),
  err_dist      = err_dists,
  contamination = contamination_levels,
  heterosked    = heterosked_types,
  stringsAsFactors = FALSE
) |>
  mutate(
    scenario_id  = row_number(),
    beta_true    = map(beta_name, ~ betas[[.x]]),
    contam_scale = 10
  )

# -----------------------------------------------------------------------------
# Build clean experiment grid
# -----------------------------------------------------------------------------

# Pinball: no hyperparameters
grid_pinball <- scenario_grid |>
  tidyr::crossing(
    tau  = taus,
    loss = losses_pinball
  ) |>
  mutate(
    k_param = NA_real_,
    c_param = NA_real_
  )

# QHuber: only k values, no c
grid_qhuber <- scenario_grid |>
  tidyr::crossing(
    tau     = taus,
    loss    = losses_qhuber,
    k_param = qhuber_k_values
  ) |>
  mutate(
    c_param = NA_real_
  )

# LogCosh: only c values, no k
grid_logcosh <- scenario_grid |>
  tidyr::crossing(
    tau     = taus,
    loss    = losses_logcosh,
    c_param = logcosh_c_values
  ) |>
  mutate(
    k_param = NA_real_
  )

# Combine all grids
experiment_grid <- bind_rows(
  grid_pinball,
  grid_qhuber,
  grid_logcosh
) |>
  mutate(
    experiment_id = row_number()
  )

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

cat("\n=== PILOT SIMULATION CONFIGURATION ===\n")
cat("Scenarios:              ", nrow(scenario_grid), "\n")
cat("  - Sample sizes:       ", paste(ns, collapse = ", "), "\n")
cat("  - Error dists:        ", paste(err_dists, collapse = ", "), "\n")
cat("  - Contamination:      ", paste(contamination_levels, collapse = ", "), "\n")
cat("  - Heterosked:         ", paste(heterosked_types, collapse = ", "), "\n")
cat("Tau levels:             ", paste(taus, collapse = ", "), "\n")
cat("Total experiments:      ", nrow(experiment_grid), "\n")
cat("Replications/experiment:", n_rep, "\n")
cat("Total fits:             ", nrow(experiment_grid) * n_rep, "\n")
cat("Estimated time:         ~1-2 minutes\n")
cat("======================================\n\n")

# Print experiment breakdown
cat("Experiment breakdown:\n")
cat("  Pinball:  ", sum(experiment_grid$loss == "pinball"), "experiments\n")
cat("  QHuber:   ", sum(experiment_grid$loss == "qhuber"), "experiments\n")
cat("  LogCosh:  ", sum(experiment_grid$loss == "logcosh"), "experiments\n")
cat("\n")

# Show first few experiments
cat("First few experiments:\n")
print(head(experiment_grid %>% 
             select(experiment_id, scenario_id, err_dist, contamination, 
                    heterosked, tau, loss, k_param, c_param), 10))

# -----------------------------------------------------------------------------
# Save configuration
# -----------------------------------------------------------------------------

save(
  scenario_grid, 
  experiment_grid, 
  n_rep,
  betas,
  ns,
  err_dists,
  taus,
  contamination_levels,
  contam_scales,
  heterosked_types,
  file = "simulation_config_SMALL.RData"
)

cat("\nSaved PILOT config to 'simulation_config_SMALL.RData'\n")