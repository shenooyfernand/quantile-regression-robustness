# =============================================================================
# Summarize Simulation Results with Monte Carlo Standard Errors
# Following Morris et al. (2019) Section 5.2 and Table 6
# =============================================================================

summarize_results <- function(results) {
  
  summary <- results %>%
    dplyr::group_by(
      experiment_id, scenario_id,
      n, beta_name, err_dist, 
      contamination, contam_scale, heterosked,
      tau, loss, k_param, c_param
    ) %>%
    dplyr::summarise(
      # Count valid replications
      n_reps = sum(!is.na(bias)),
      
      # ========================================================================
      # PERFORMANCE MEASURES (Averaged Across Replications)
      # ========================================================================
      mean_bias  = mean(bias, na.rm = TRUE),
      mean_mse   = mean(mse, na.rm = TRUE),
      mean_mae   = mean(mae, na.rm = TRUE),
      mean_time  = mean(computation_time, na.rm = TRUE),
      conv_rate  = mean(converged, na.rm = TRUE),
      mean_iter  = mean(iterations, na.rm = TRUE),
      
      # ========================================================================
      # MONTE CARLO STANDARD ERRORS
      # Following Morris et al. (2019) Table 6, page 13
      # ========================================================================
      
      # MC SE for bias: SD(bias_i) / sqrt(n_reps)
      # Formula: sqrt(Var(theta_hat)) / sqrt(n_sim)
      mcse_bias = sd(bias, na.rm = TRUE) / sqrt(sum(!is.na(bias))),
      
      # MC SE for MSE: 
      # Formula: sqrt(sum((mse_i - mean_mse)^2) / (n_sim * (n_sim - 1)))
      mcse_mse = {
        valid_mse <- mse[!is.na(mse)]
        n_valid <- length(valid_mse)
        if (n_valid > 1) {
          mean_mse_val <- mean(valid_mse)
          sqrt(sum((valid_mse - mean_mse_val)^2) / (n_valid * (n_valid - 1)))
        } else {
          NA_real_
        }
      },
      
      # MC SE for MAE: SD(mae_i) / sqrt(n_reps)
      mcse_mae = sd(mae, na.rm = TRUE) / sqrt(sum(!is.na(mae))),
      
      # MC SE for convergence rate (proportion):
      # Formula: sqrt(p * (1-p) / n_sim)
      mcse_conv = {
        p <- mean(converged, na.rm = TRUE)
        n_valid <- sum(!is.na(converged))
        if (n_valid > 0 && !is.na(p)) {
          sqrt(p * (1 - p) / n_valid)
        } else {
          NA_real_
        }
      },
      
      # MC SE for computation time
      mcse_time = sd(computation_time, na.rm = TRUE) / sqrt(sum(!is.na(computation_time))),
      
      # MC SE for iterations
      mcse_iter = sd(iterations, na.rm = TRUE) / sqrt(sum(!is.na(iterations))),
      
      .groups = "drop"
    )
  
  # Add quality indicators for Monte Carlo precision
  summary <- summary %>%
    mutate(
      # Check if MC SEs meet common targets
      bias_precision = case_when(
        is.na(mcse_bias) ~ "Unknown",
        mcse_bias < 0.005 ~ "Excellent",
        mcse_bias < 0.01 ~ "Good",
        mcse_bias < 0.02 ~ "Acceptable",
        TRUE ~ "Increase n_rep"
      ),
      
      conv_precision = case_when(
        is.na(mcse_conv) ~ "Unknown",
        mcse_conv < 0.005 ~ "Excellent",
        mcse_conv < 0.01 ~ "Good",
        mcse_conv < 0.02 ~ "Acceptable",
        TRUE ~ "Increase n_rep"
      ),
      
      # Overall assessment
      overall_precision = case_when(
        bias_precision %in% c("Excellent", "Good") & 
          conv_precision %in% c("Excellent", "Good") ~ "Adequate",
        TRUE ~ "Review"
      )
    )
  
  return(summary)
}


# =============================================================================
# Diagnostic Function: Check Monte Carlo SE Adequacy
# Following Morris et al. (2019) Section 5.3
# =============================================================================

check_monte_carlo_precision <- function(summary_results, 
                                        target_mcse_bias = 0.01,
                                        target_mcse_conv = 0.01) {
  
  cat("\n")
  cat(rep("=", 75), "\n", sep = "")
  cat("       MONTE CARLO STANDARD ERROR DIAGNOSTIC\n")
  cat("       Following Morris et al. (2019) Section 5.3\n")
  cat(rep("=", 75), "\n\n", sep = "")
  
  # Get n_rep (should be same across all experiments)
  n_rep_values <- unique(summary_results$n_reps)
  if (length(n_rep_values) == 1) {
    cat("Current simulation: n_rep =", n_rep_values, "replications per experiment\n\n")
  } else {
    cat("Warning: Different n_rep values detected:", 
        paste(n_rep_values, collapse = ", "), "\n\n")
  }
  
  # -------------------------------------------------------------------------
  # Check BIAS Monte Carlo SE
  # -------------------------------------------------------------------------
  cat("BIAS PRECISION:\n")
  max_mcse_bias <- max(abs(summary_results$mcse_bias), na.rm = TRUE)
  mean_mcse_bias <- mean(abs(summary_results$mcse_bias), na.rm = TRUE)
  
  cat("  Maximum MC SE:  ", sprintf("%.5f", max_mcse_bias), "\n")
  cat("  Mean MC SE:     ", sprintf("%.5f", mean_mcse_bias), "\n")
  cat("  Target:          < ", target_mcse_bias, "\n")
  cat("  Status:         ", 
      ifelse(max_mcse_bias < target_mcse_bias, 
             "✓ ADEQUATE", 
             "✗ INCREASE n_rep"), 
      "\n\n")
  
  # -------------------------------------------------------------------------
  # Check CONVERGENCE RATE Monte Carlo SE
  # -------------------------------------------------------------------------
  cat("CONVERGENCE RATE PRECISION:\n")
  max_mcse_conv <- max(summary_results$mcse_conv, na.rm = TRUE)
  mean_mcse_conv <- mean(summary_results$mcse_conv, na.rm = TRUE)
  
  cat("  Maximum MC SE:  ", sprintf("%.5f (%.2f%%)", max_mcse_conv, max_mcse_conv*100), "\n")
  cat("  Mean MC SE:     ", sprintf("%.5f (%.2f%%)", mean_mcse_conv, mean_mcse_conv*100), "\n")
  cat("  Target:          < ", target_mcse_conv, sprintf("(%.1f%%)", target_mcse_conv*100), "\n")
  cat("  Status:         ", 
      ifelse(max_mcse_conv < target_mcse_conv, 
             "✓ ADEQUATE", 
             "✗ INCREASE n_rep"), 
      "\n\n")
  
  # -------------------------------------------------------------------------
  # Check MSE Monte Carlo SE
  # -------------------------------------------------------------------------
  cat("MSE PRECISION:\n")
  max_mcse_mse <- max(summary_results$mcse_mse, na.rm = TRUE)
  mean_mcse_mse <- mean(summary_results$mcse_mse, na.rm = TRUE)
  
  cat("  Maximum MC SE:  ", sprintf("%.5f", max_mcse_mse), "\n")
  cat("  Mean MC SE:     ", sprintf("%.5f", mean_mcse_mse), "\n\n")
  
  # -------------------------------------------------------------------------
  # Recommendations if inadequate
  # -------------------------------------------------------------------------
  needs_increase <- (max_mcse_bias >= target_mcse_bias) || 
    (max_mcse_conv >= target_mcse_conv)
  
  if (needs_increase) {
    current_n <- n_rep_values[1]
    
    cat(rep("-", 75), "\n", sep = "")
    cat("RECOMMENDATIONS:\n\n")
    
    if (max_mcse_bias >= target_mcse_bias) {
      # Estimate SD(bias) from current MC SE
      sd_bias_est <- max_mcse_bias * sqrt(current_n)
      required_n_bias <- ceiling((sd_bias_est / target_mcse_bias)^2)
      
      cat("For bias precision:\n")
      cat("  Current n_rep:  ", current_n, "\n")
      cat("  Required n_rep: ", required_n_bias, "\n")
      cat("  Increase factor:", sprintf("%.1fx", required_n_bias / current_n), "\n\n")
    }
    
    if (max_mcse_conv >= target_mcse_conv) {
      # Use worst case: p = 0.5 for conservative estimate
      required_n_conv <- ceiling(0.25 / target_mcse_conv^2)
      
      cat("For convergence rate precision:\n")
      cat("  Current n_rep:  ", current_n, "\n")
      cat("  Required n_rep: ", required_n_conv, "(worst case)\n")
      cat("  Increase factor:", sprintf("%.1fx", required_n_conv / current_n), "\n\n")
    }
    
    cat("Action: Re-run simulation with increased n_rep\n")
    cat(rep("-", 75), "\n", sep = "")
    
  } else {
    cat(rep("-", 75), "\n", sep = "")
    cat("✓ SIMULATION PRECISION IS ADEQUATE\n")
    cat("  All Monte Carlo SEs meet target precision.\n")
    cat("  Results can be reported with confidence.\n")
    cat(rep("-", 75), "\n", sep = "")
  }
  
  cat("\n")
  cat(rep("=", 75), "\n\n", sep = "")
  
  invisible(summary_results)
}


# =============================================================================
# Function to Create Publication-Ready Tables with MC SEs
# =============================================================================

format_results_table <- function(summary_results, 
                                 selected_scenarios = NULL,
                                 digits_bias = 4,
                                 digits_mse = 4,
                                 digits_conv = 3) {
  
  # Filter scenarios if specified
  if (!is.null(selected_scenarios)) {
    summary_results <- summary_results %>%
      filter(scenario_id %in% selected_scenarios)
  }
  
  # Create formatted display columns
  formatted <- summary_results %>%
    mutate(
      # Format: "estimate (MC SE)"
      Bias = sprintf(paste0("%.", digits_bias, "f (%.", digits_bias, "f)"), 
                     mean_bias, mcse_bias),
      MSE  = sprintf(paste0("%.", digits_mse, "f (%.", digits_mse, "f)"), 
                     mean_mse, mcse_mse),
      MAE  = sprintf(paste0("%.", digits_mse, "f (%.", digits_mse, "f)"), 
                     mean_mae, mcse_mae),
      Conv = sprintf(paste0("%.", digits_conv, "f (%.", digits_conv, "f)"), 
                     conv_rate, mcse_conv),
      Time = sprintf("%.2f (%.2f)", mean_time, mcse_time)
    ) %>%
    select(
      experiment_id, scenario_id, 
      n, beta_name, err_dist, contamination, heterosked, tau, 
      loss, k_param, c_param,
      Bias, MSE, MAE, Conv, Time
    ) %>%
    arrange(scenario_id, tau, loss)
  
  return(formatted)
}