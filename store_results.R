store_results <- function(results, file = "simulation_results_raw.RData") {
  save(results, file = file)
  invisible(results)
}