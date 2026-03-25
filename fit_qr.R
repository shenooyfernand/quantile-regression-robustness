fit_qr <- function(X, y, tau, loss, k_param = NULL, c_param = NULL) {
  # X: n x p
  # y: length n
  # tau: quantile
  # loss: "pinball", "qhuber", "logcosh"
  # k_param: k for qhuber
  # c_param: c for logcosh
  
  n <- nrow(X)
  p <- ncol(X)
  
  obj <- function(beta) {
    r <- as.numeric(y - X %*% beta)
    ld <- loss_dispatch(loss, r, X, tau, k = k_param, c = c_param)
    sum(ld$loss)
  }
  
  grad <- function(beta) {
    r <- as.numeric(y - X %*% beta)
    ld <- loss_dispatch(loss, r, X, tau, k = k_param, c = c_param)
    ld$grad
  }
  
  beta0 <- rep(0, p)
  
  out <- tryCatch(
    optim(
      par     = beta0,
      fn      = obj,
      gr      = grad,
      method  = "L-BFGS-B",
      control = list(maxit = 1000)
    ),
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(list(
      beta       = rep(NA_real_, p),
      converged  = FALSE,
      iterations = NA_integer_
    ))
  }
  
  beta_hat <- out$par
  
  # If any NA in coefficients, mark as non-converged
  if (any(is.na(beta_hat))) {
    return(list(
      beta       = rep(NA_real_, p),
      converged  = FALSE,
      iterations = NA_integer_
    ))
  }
  
  list(
    beta       = as.numeric(beta_hat),
    converged  = (out$convergence == 0),
    iterations = out$counts[["function"]]
  )
}