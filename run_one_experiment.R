
run_one_experiment <- function(exp_row, n_rep) {
  purrr::map_dfr(
    seq_len(n_rep),
    ~ run_one_rep(exp_row, rep_id = .x)
  )
}