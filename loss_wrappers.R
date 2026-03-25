# =============================================================================
# Loss Function Wrappers for Residual-Based Computation
# =============================================================================
# These wrappers accept pre-computed residuals (r) instead of beta, x, y
# This matches the interface expected by loss_dispatch()

# Pinball Loss Wrappers
# =============================================================================
pinball_loss <- function(r, tau) {
  # r: residuals (y - X %*% beta), already computed
  # tau: quantile level
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must be between 0 and 1")
  }
  
  # Pinball loss: rho_tau(r) = r * (tau - I(r < 0))
  loss <- r * (tau - (r < 0))
  loss
}

pinball_grad <- function(r, X, tau) {
  # Gradient of pinball loss
  # psi(r) = tau if r >= 0, else (tau - 1)
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must be between 0 and 1")
  }
  
  psi <- ifelse(r >= 0, tau, tau - 1)
  
  # gradient wrt beta: -X^T * psi
  grad <- -crossprod(X, psi)
  as.numeric(grad)
}


# =============================================================================
# Quantile Huber Loss Wrappers (Koenker-consistent)
# Residual convention: r = y - X %*% beta
# Targets: P(r <= 0) = tau
# =============================================================================

qhuber_loss <- function(r, tau, k) {
  if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
  if (k <= 0) stop("k must be positive")
  
  # Thresholds (Koenker convention)
  l <- (tau - 1) * k   # negative bound
  u <- tau * k         # positive bound
  
  loss <- numeric(length(r))
  
  low  <- r < l
  mid  <- r >= l & r <= u
  high <- r > u
  
  # Lower tail: linear with slope (tau - 1)
  # Since r<0 here, this is positive loss: (tau-1)*r - 0.5*k*(tau-1)^2
  loss[low] <- (1-tau) * abs(r[low]) - 0.5 * k * (1-tau)^2
  
  # Middle: quadratic
  loss[mid] <- 0.5 * r[mid]^2 / k
  
  # Upper tail: linear with slope tau
  loss[high] <- tau * r[high] - 0.5 * k * tau^2
  
  loss
}

qhuber_grad <- function(r, X, tau, k) {
  if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
  if (k <= 0) stop("k must be positive")
  
  l <- (tau - 1) * k
  u <- tau * k
  
  psi <- numeric(length(r))
  
  low  <- r < l
  mid  <- r >= l & r <= u
  high <- r > u
  
  # Score / derivative w.r.t r
  psi[low]  <- (tau - 1)
  psi[mid]  <- r[mid] / k
  psi[high] <- tau
  
  # r = y - X beta  => dr/dbeta = -X, so grad = -X^T psi
  grad <- -crossprod(X, psi)
  as.numeric(grad)
}


# LogCosh Quantile Loss Wrappers
# =============================================================================
logcosh_loss <- function(r, tau, c) {
  if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
  if (c <= 0) stop("c must be positive")
  
  x <- c * r
  
  # Numerically stable log(cosh(x))
  # For |x| > 20: log(cosh(x)) b	 |x| - log(2)
  abs_x <- abs(x)
  
  log_cosh <- ifelse(
    abs_x > 20,
    abs_x - log(2),  # Asymptotic approximation
    log(cosh(x))     # Direct computation
  )
  
  loss <- (1/(2*c)) * log_cosh + (tau - 0.5) * r
  loss
}

logcosh_grad <- function(r, X, tau, c) {
  if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
  if (c <= 0) stop("c must be positive")
  
  x <- c * r
  
  # Numerically stable tanh(x)
  # For |x| > 20: tanh(x) b	 sign(x)
  tanh_x <- ifelse(
    abs(x) > 20,
    sign(x),   # Asymptotic value
    tanh(x)    # Direct computation
  )
  
  psi <- 0.5 * tanh_x + (tau - 0.5)
  
  grad <- -crossprod(X, psi)
  as.numeric(grad)
}