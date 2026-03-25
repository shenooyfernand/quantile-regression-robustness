run_one_rep <- function(exp_row, rep_id) {
  n           <- exp_row$n
  beta_true   <- as.numeric(exp_row$beta_true[[1]])
  err_dist    <- exp_row$err_dist
  tau         <- exp_row$tau
  loss        <- exp_row$loss
  k_param     <- exp_row$k_param
  c_param     <- exp_row$c_param
  contamination <- exp_row$contamination
  contam_scale  <- exp_row$contam_scale
  heterosked    <- exp_row$heterosked
  
  # 1) simulate data
  dat <- simulate_data(
    n             = n, 
    beta_true     = beta_true, 
    err_dist      = err_dist,
    tau           = tau,
    contamination = contamination,
    contam_scale  = contam_scale,
    heterosked    = heterosked
  )
  X <- dat$X
  y <- dat$y
  
  # 2) fit model
  t_start <- Sys.time()
  fit <- fit_qr(
    X       = X,
    y       = y,
    tau     = tau,
    loss    = loss,
    k_param = if (!is.na(k_param)) k_param else NULL,
    c_param = if (!is.na(c_param)) c_param else NULL
  )
  t_end <- Sys.time()
  
  beta_hat   <- fit$beta
  converged  <- isTRUE(fit$converged)
  iterations <- if (!is.null(fit$iterations)) fit$iterations else NA_integer_
  
  # 3) metrics
  bias        <- compute_bias(beta_hat, beta_true)
  mse         <- compute_mse(beta_hat, beta_true)
  mae         <- compute_mae(beta_hat, beta_true)
  comp_time   <- compute_computation_time(t_start, t_end)
  
  # REMOVED: coverage_95 and ci_width (they would be NA anyway)
  
  tibble::tibble(
    experiment_id     = exp_row$experiment_id,
    scenario_id       = exp_row$scenario_id,
    rep               = rep_id,
    n                 = n,
    beta_name         = exp_row$beta_name,
    err_dist          = err_dist,
    contamination     = contamination,
    contam_scale      = contam_scale,
    heterosked        = heterosked,
    tau               = tau,
    loss              = loss,
    k_param           = k_param,
    c_param           = c_param,
    bias              = bias,
    mse               = mse,
    mae               = mae,
    computation_time  = comp_time,
    converged         = converged,
    iterations        = iterations
  )
}