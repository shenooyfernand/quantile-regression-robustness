## Robustness in Quantile Regression: An Influence Function 
**Tools:** R

### Abstract
Quantile regression models different parts of the conditional 
distribution of a response variable and is widely used to study 
tail behaviour and heteroskedasticity. This thesis investigates 
how the choice of loss function affects robustness and estimation 
accuracy in quantile regression.

### Loss Functions Compared
- Pinball loss (classical standard)
- Quantile Huber loss
- Logcosh-based smooth loss

### Methods
- Influence Function theory to analytically assess robustness
  of each loss function
- Large-scale Monte Carlo simulation study:
  - 1,000 replications per scenario
  - Varying sample sizes and quantile levels (0.25, 0.50, 0.75)
  - Multiple contamination proportions
  - Heteroskedastic variance structures
  - Multiple error distributions (normal, heavy-tailed)
- Performance evaluated using MSE of regression coefficients, 
  bias, variability, and computational measures

### Key Findings
- Quantile Huber loss consistently produced the lowest MSE 
  across simulation settings
- Robust loss functions outperform pinball loss under response 
  contamination and heavy-tailed error distributions
- Influence function analysis explains this: the piecewise 
  structure of the Huber score function controls extreme 
  residuals without distorting estimation near the centre
- Quantile Huber loss is the preferred choice when contamination 
  or heavy-tailed errors are expected, provided the tuning 
  parameter is selected carefully

### Repository Contents
- R simulation scripts
- Simulation output data
