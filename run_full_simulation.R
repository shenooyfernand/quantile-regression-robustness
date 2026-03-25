run_full_simulation <- function(experiment_grid, n_rep,
                                parallel = TRUE,
                                seed = TRUE) {
  
  if (parallel) {
    results <- furrr::future_map_dfr(
      1:nrow(experiment_grid),
      function(i) {
        exp_row <- experiment_grid[i, ]
        run_one_experiment(exp_row, n_rep = n_rep)
      },
      .progress = TRUE,
      .options  = furrr::furrr_options(seed = seed)
    )
  } else {
    results <- purrr::map_dfr(
      1:nrow(experiment_grid),
      function(i) {
        exp_row <- experiment_grid[i, ]
        run_one_experiment(exp_row, n_rep = n_rep)
      }
    )
  }
  
  results
}
