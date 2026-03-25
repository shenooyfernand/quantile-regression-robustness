save_simulation_results <- function(results, summary_results, 
                                    experiment_grid = NULL,
                                    prefix = "simulation",
                                    format = c("excel", "csv", "both")) {
  
  format <- match.arg(format)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create output directory
  dir.create("simulation_output", showWarnings = FALSE)
  
  # ========================================================================
  # Summary with MC SEs (include all new columns)
  # ========================================================================
  summary_selected <- summary_results %>%
    dplyr::select(
      experiment_id, scenario_id, n, beta_name, err_dist,
      contamination, contam_scale, heterosked,
      tau, loss, k_param, c_param,
      n_reps,  # NEW: number of replications
      mean_bias, mcse_bias,  # NEW: MC SE for bias
      mean_mse, mcse_mse,    # NEW: MC SE for MSE
      mean_mae, mcse_mae,    # NEW: MC SE for MAE
      mean_time, mcse_time,  # NEW: MC SE for time
      conv_rate, mcse_conv,  # NEW: MC SE for convergence
      mean_iter, mcse_iter,  # NEW: MC SE for iterations
      bias_precision, conv_precision, overall_precision  # NEW: quality flags
    )
  
  # Best performers by MSE for each scenario/tau
  best_by_scenario <- summary_results %>%
    dplyr::group_by(scenario_id, tau) %>%
    dplyr::slice_min(mean_mse, n = 3, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(scenario_id, tau, mean_mse)
  
  # Performance by loss type
  perf_by_loss <- summary_results %>%
    dplyr::group_by(loss, k_param, c_param) %>%
    dplyr::summarise(
      n_experiments = n(),
      avg_mse = mean(mean_mse, na.rm = TRUE),
      avg_mae = mean(mean_mae, na.rm = TRUE),
      avg_bias = mean(abs(mean_bias), na.rm = TRUE),  # Absolute bias
      avg_conv_rate = mean(conv_rate, na.rm = TRUE),
      # NEW: Average MC SEs
      avg_mcse_bias = mean(mcse_bias, na.rm = TRUE),
      avg_mcse_mse = mean(mcse_mse, na.rm = TRUE),
      avg_mcse_conv = mean(mcse_conv, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(loss, k_param, c_param)
  
  # Performance under contamination
  perf_by_contam <- summary_results %>%
    dplyr::group_by(contamination, loss) %>%
    dplyr::summarise(
      avg_mse = mean(mean_mse, na.rm = TRUE),
      avg_mae = mean(mean_mae, na.rm = TRUE),
      avg_bias = mean(abs(mean_bias), na.rm = TRUE),
      n_experiments = n(),
      # NEW: Max MC SE (worst case precision)
      max_mcse_bias = max(mcse_bias, na.rm = TRUE),
      max_mcse_mse = max(mcse_mse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(contamination, avg_mse)
  
  # Performance under heteroskedasticity
  perf_by_hetero <- summary_results %>%
    dplyr::group_by(heterosked, loss) %>%
    dplyr::summarise(
      avg_mse = mean(mean_mse, na.rm = TRUE),
      avg_mae = mean(mean_mae, na.rm = TRUE),
      avg_bias = mean(abs(mean_bias), na.rm = TRUE),
      n_experiments = n(),
      # NEW: Max MC SE
      max_mcse_bias = max(mcse_bias, na.rm = TRUE),
      max_mcse_mse = max(mcse_mse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(heterosked, avg_mse)
  
  # ========================================================================
  # Monte Carlo Precision Summary
  # ========================================================================
  mc_precision_summary <- summary_results %>%
    dplyr::summarise(
      n_rep = unique(n_reps)[1],
      max_mcse_bias = max(mcse_bias, na.rm = TRUE),
      mean_mcse_bias = mean(mcse_bias, na.rm = TRUE),
      max_mcse_mse = max(mcse_mse, na.rm = TRUE),
      mean_mcse_mse = mean(mcse_mse, na.rm = TRUE),
      max_mcse_conv = max(mcse_conv, na.rm = TRUE),
      mean_mcse_conv = mean(mcse_conv, na.rm = TRUE),
      pct_adequate_bias = mean(bias_precision %in% c("Excellent", "Good")) * 100,
      pct_adequate_conv = mean(conv_precision %in% c("Excellent", "Good")) * 100
    )
  
  # ========================================================================
  # Save to Excel
  # ========================================================================
  if (format %in% c("excel", "both")) {
    library(writexl)
    
    excel_list <- list(
      README = data.frame(
        Info = c(
          "Monte Carlo Standard Errors (MC SE) are included",
          "MC SEs quantify simulation uncertainty",
          "Following Morris et al. (2019) Statistics in Medicine",
          "Columns with 'mcse_' prefix show Monte Carlo SEs",
          "See MC_Precision sheet for adequacy assessment",
          paste("Simulation run:", timestamp),
          paste("n_rep =", mc_precision_summary$n_rep)
        )
      ),
      MC_Precision = mc_precision_summary,
      Summary_with_MCSE = summary_selected,
      Best_By_Scenario = best_by_scenario,
      Perf_By_Loss = perf_by_loss,
      Perf_By_Contamination = perf_by_contam,
      Perf_By_Heterosked = perf_by_hetero,
      Summary_Full = summary_results,
      Raw_Sample = head(results, 2000)
    )
    
    if (!is.null(experiment_grid)) {
      excel_list$Config = experiment_grid
    }
    
    excel_file <- paste0("simulation_output/", prefix, "_results_", timestamp, ".xlsx")
    write_xlsx(excel_list, path = excel_file)
    cat("Excel results saved to:", excel_file, "\n")
  }
  
  # ========================================================================
  # Save to CSV
  # ========================================================================
  if (format %in% c("csv", "both")) {
    csv_prefix <- paste0("simulation_output/", prefix, "_", timestamp)
    
    readr::write_csv(mc_precision_summary, paste0(csv_prefix, "_mc_precision.csv"))
    readr::write_csv(summary_selected, paste0(csv_prefix, "_summary.csv"))
    readr::write_csv(best_by_scenario, paste0(csv_prefix, "_best.csv"))
    readr::write_csv(perf_by_loss, paste0(csv_prefix, "_by_loss.csv"))
    readr::write_csv(perf_by_contam, paste0(csv_prefix, "_by_contamination.csv"))
    readr::write_csv(perf_by_hetero, paste0(csv_prefix, "_by_heterosked.csv"))
    readr::write_csv(results, paste0(csv_prefix, "_raw.csv"))
    
    cat("CSV results saved with prefix:", csv_prefix, "\n")
  }
  
  # ========================================================================
  # Save RData
  # ========================================================================
  rdata_file <- paste0("simulation_output/", prefix, "_", timestamp, ".RData")
  save(results, summary_results, mc_precision_summary, file = rdata_file)
  cat("RData saved to:", rdata_file, "\n")
  
  invisible(list(
    summary = summary_selected,
    best = best_by_scenario,
    mc_precision = mc_precision_summary,
    timestamp = timestamp
  ))
}