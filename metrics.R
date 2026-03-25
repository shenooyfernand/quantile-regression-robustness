# =============================================================================
# Performance Metrics for Simulation Study
# =============================================================================
# These compute per-replication metrics from estimates and true values
# Monte Carlo SEs are computed later in summarize_results.R

compute_bias <- function(beta_hat, beta_true) {
  if (any(is.na(beta_hat))) return(NA_real_)
  mean(beta_hat - beta_true)
}

compute_mse <- function(beta_hat, beta_true) {
  if (any(is.na(beta_hat))) return(NA_real_)
  diff <- beta_hat - beta_true
  mean(diff^2)
}

compute_mae <- function(beta_hat, beta_true) {
  if (any(is.na(beta_hat))) return(NA_real_)
  diff <- beta_hat - beta_true
  mean(abs(diff))
}

# These functions are not currently used since se_beta is always NULL
# Keeping them for potential future use with real data analysis
compute_coverage_95 <- function(beta_hat, beta_true, se_beta = NULL) {
  if (is.null(se_beta) || any(is.na(beta_hat))) return(NA_real_)
  z <- 1.96
  lower <- beta_hat - z * se_beta
  upper <- beta_hat + z * se_beta
  covered <- (beta_true >= lower) & (beta_true <= upper)
  mean(covered)
}

compute_ci_width <- function(se_beta = NULL) {
  if (is.null(se_beta)) return(NA_real_)
  z <- 1.96
  mean(2 * z * se_beta)
}

compute_convergence_stats <- function(converged, iterations = NA_integer_) {
  list(
    convergence = as.logical(converged),
    iterations  = iterations
  )
}

compute_computation_time <- function(t_start, t_end) {
  as.numeric(difftime(t_end, t_start, units = "secs"))
}