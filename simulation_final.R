# =============================================================================
# Main Simulation Runner with Monte Carlo SE Diagnostics
# =============================================================================

# Load all files
source("loss_wrappers.R")
source("lossdispather.R")
source("fit_qr.R")
source("simulate_data.R")
source("one_rep.R")
source("run_one_experiment.R")
source("run_full_simulation.R")
source("metrics.R")
source("summarize_results.R")
source("store_results.R")
source("save_results.R")

# Load config
load("simulation_config.RData")

# Print configuration info
cat("\n=== Starting Full Simulation ===\n")
cat("Total experiments:", nrow(experiment_grid), "\n")
cat("Replications:", n_rep, "\n")
cat("Total fits:", nrow(experiment_grid) * n_rep, "\n")
cat("\nMonte Carlo Precision Targets:\n")
cat("  MCSE(bias):       <", config_info$mc_precision$target_mcse_bias, "\n")
cat("  MCSE(conv):       <", config_info$mc_precision$target_mcse_conv, "\n")
cat("  Expected MCSE(bias):", sprintf("%.4f", config_info$mc_precision$achieved_mcse_bias), "\n")
cat("  Expected MCSE(conv):", sprintf("%.4f", config_info$mc_precision$achieved_mcse_conv), "\n")
cat("\n")

# Run simulation
results <- run_full_simulation(experiment_grid, n_rep = n_rep, parallel = TRUE)

# Summarize with Monte Carlo SEs
cat("\n=== Summarizing Results (with Monte Carlo SEs) ===\n")
summary_results <- summarize_results(results)

# ============================================================================
# Check Monte Carlo Precision
# ============================================================================
cat("\n=== Checking Monte Carlo Precision ===\n")
check_monte_carlo_precision(
  summary_results,
  target_mcse_bias = config_info$mc_precision$target_mcse_bias,
  target_mcse_conv = config_info$mc_precision$target_mcse_conv
)

# Quick preview of results
cat("\n=== Quick Preview (First 10 Experiments) ===\n")
preview <- summary_results %>%
  head(10) %>%
  select(experiment_id, loss, k_param, c_param, tau,
         mean_bias, mcse_bias, mean_mse, mcse_mse, 
         conv_rate, mcse_conv)

print(preview, n = 10)

# Save all results
cat("\n=== Saving Results ===\n")
save_simulation_results(
  results = results,
  summary_results = summary_results,
  experiment_grid = experiment_grid,
  prefix = "FULL",
  format = "both"
)

cat("\n=== Full Simulation Complete ===\n")
cat("Check 'simulation_output/' folder for results\n")
cat("Key files:\n")
cat("  - *_mc_precision.csv : Monte Carlo SE adequacy report\n")
cat("  - *_summary.csv      : Summary with MC SEs\n")
cat("  - *.xlsx             : Excel workbook with all results\n\n")

