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
load("simulation_config_SMALL.RData")

# Run simulation
cat("\n=== Running Pilot Simulation ===\n")
results <- run_full_simulation(experiment_grid, n_rep = n_rep, parallel = TRUE)

# Summarize
cat("\n=== Summarizing Results ===\n")
summary_results <- summarize_results(results)

# Quick view
cat("\n=== Quick Summary ===\n")
summary_results %>%
  dplyr::select(experiment_id, loss, k_param, c_param,
                mean_bias, mean_mse, mean_mae, conv_rate) %>%
  print()

# Save all results
cat("\n=== Saving Results ===\n")
save_simulation_results(
  results = results,
  summary_results = summary_results,
  experiment_grid = experiment_grid,
  prefix = "PILOT",
  format = "both"
)

cat("\n=== Pilot Simulation Complete ===\n")

