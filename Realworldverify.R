rm(list=ls())
# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------
library(glmnet)
library(dplyr)

# ------------------------------------------------------------
# 2. Load training data
# ------------------------------------------------------------
train <- read.csv("C:/Users/User/Downloads/house-prices-advanced-regression-techniques/train.csv")

# Remove ID
train$Id <- NULL

# ------------------------------------------------------------
# 3. Convert character variables to factor
# ------------------------------------------------------------
char_vars <- names(train)[sapply(train, is.character)]
for (v in char_vars) {
  train[[v]] <- as.factor(train[[v]])
}

# ------------------------------------------------------------
# 4. Handle missing values
# ------------------------------------------------------------

# ---- 4a. Replace NA in factor variables with "None"
factor_vars <- names(train)[sapply(train, is.factor)]

for (v in factor_vars) {
  if (any(is.na(train[[v]]))) {
    train[[v]] <- addNA(train[[v]])
    levels(train[[v]])[is.na(levels(train[[v]]))] <- "None"
  }
}

# ---- 4b. Median imputation for numeric variables
numeric_vars <- names(train)[sapply(train, is.numeric)]

for (v in numeric_vars) {
  if (any(is.na(train[[v]]))) {
    train[[v]][is.na(train[[v]])] <- median(train[[v]], na.rm = TRUE)
  }
}

# ------------------------------------------------------------
# 5. Final check: there must be NO NA remaining
# ------------------------------------------------------------
if (sum(is.na(train)) > 0) {
  stop("There are still missing values!")
}

# ------------------------------------------------------------
# 6. Define response and predictors
# ------------------------------------------------------------
y <- log(train$SalePrice)
X_raw <- train %>% select(-SalePrice)

# ------------------------------------------------------------
# 7. Create model matrix (dummy encoding)
# ------------------------------------------------------------
X <- model.matrix(~ ., data = X_raw)[, -1]  # remove intercept

# Confirm dimensions match
cat("Rows in X:", nrow(X), "\n")
cat("Length of y:", length(y), "\n")

# ------------------------------------------------------------
# 8. Fit cross-validated LASSO
# ------------------------------------------------------------
set.seed(123)

cv_lasso <- cv.glmnet(
  X,
  y,
  alpha = 1,        # LASSO penalty
  nfolds = 10
)

# Use lambda.1se for more stable model
lambda_opt <- cv_lasso$lambda.1se

# ------------------------------------------------------------
# 9. Fit final LASSO model
# ------------------------------------------------------------
lasso_model <- glmnet(X, y, alpha = 1, lambda = lambda_opt)

# Extract coefficients
coef_lasso <- coef(lasso_model)

# Get non-zero variables
selected_vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
selected_vars <- selected_vars[selected_vars != "(Intercept)"]

cat("Selected variables:\n")
print(selected_vars)


# Response
y_real <- log(train$SalePrice)

# Predictors (only selected variables)
X_selected <- X[, selected_vars]

# Add intercept column (if your fit_qr expects it)
X_selected <- cbind(Intercept = 1, X_selected)

dim(X_selected)
length(y_real)

taus <- c(0.1, 0.5, 0.9)
source("fit_qr.R")

# Standardize predictors
X_scaled <- X_selected
X_scaled[, -1] <- scale(X_scaled[, -1])
X_selected <- X_scaled
results_real <- list()

for (tau in taus) {
  
  # Pinball
  fit_pin <- fit_qr(
    X = X_selected,
    y = y_real,
    tau = tau,
    loss = "pinball"
  )
  
  # Quantile Huber (choose k consistent with simulation)
  fit_qh <- fit_qr(
    X = X_selected,
    y = y_real,
    tau = tau,
    loss = "qhuber",
    k_param = 1.0   # use same k as simulation
  )
  
  # Log-cosh (choose c consistent with simulation)
  fit_lc <- fit_qr(
    X = X_selected,
    y = y_real,
    tau = tau,
    loss = "logcosh",
    c_param = 1.0   # same c as simulation
  )
  
  results_real[[paste0("tau_", tau)]] <- list(
    pinball = fit_pin,
    qhuber  = fit_qh,
    logcosh = fit_lc
  )
}

check_loss <- function(r, tau) {
  sum(r * (tau - (r < 0)))
}

performance <- list()

for (tau in taus) {
  
  perf_tau <- list()
  
  for (loss_name in c("pinball", "qhuber", "logcosh")) {
    
    beta_hat <- results_real[[paste0("tau_", tau)]][[loss_name]]$beta
    residuals <- y_real - X_selected %*% beta_hat
    
    perf_tau[[loss_name]] <- check_loss(residuals, tau)
  }
  
  performance[[paste0("tau_", tau)]] <- perf_tau
}

performance

