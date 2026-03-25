simulate_data <- function(n, beta_true, err_dist = c("normal", "exp"), tau = 0.5,
                          contamination = 0, contam_scale = 10,
                          heterosked = c("none", "linear", "quadratic")) {
  
  err_dist <- match.arg(err_dist)
  heterosked <- match.arg(heterosked)
  p <- length(beta_true)
  
  # Generate X with p-1 columns (no intercept yet)
  Sigma <- outer(1:(p-1), 1:(p-1), function(j, k) 0.5^(abs(j - k)))
  X_no_int <- MASS::mvrnorm(n = n, mu = rep(0, p-1), Sigma = Sigma)
  X <- cbind(1, X_no_int)
  
  # Generate base errors with Q_τ(ε) = 0
  if (err_dist == "normal") {
    eps_raw <- rnorm(n, mean = 0, sd = 1)
    q_tau <- qnorm(tau, mean = 0, sd = 1)
    eps <- eps_raw - q_tau
    
  } else if (err_dist == "exp") {
    Z <- rexp(n, rate = 1)
    q_tau <- qexp(tau, rate = 1)
    eps <- Z - q_tau
  }
  
  # Apply heteroskedasticity
  if (heterosked == "linear") {
    # Variance proportional to |X[,2]|
    scale_factor <- sqrt(1 + abs(X[, 2]))
    eps <- eps * scale_factor
    
  } else if (heterosked == "quadratic") {
    # Variance proportional to X[,2]^2
    scale_factor <- sqrt(1 + X[, 2]^2)
    eps <- eps * scale_factor
  }
  # heterosked == "none": no scaling
  
  # Add contamination
  if (contamination > 0) {
    n_contam <- floor(n * contamination)
    if (n_contam > 0) {
      contam_idx <- sample(n, n_contam)
      # Additive outliers
      eps[contam_idx] <- eps[contam_idx] + rnorm(n_contam, 0, contam_scale)
    }
  }
  
  # Generate response
  y <- as.numeric(X %*% beta_true + eps)
  
  list(X = X, y = y)
}