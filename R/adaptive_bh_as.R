

# Define the AS method based on the provided description
adaptive_bh_as <- function(p_values, q, delta = 0.01) {
  n <- length(p_values)

  # Step 1: Find the optimal lambda using the stopping rule
  lambda_grid <- seq(q, 1, by = delta)
  pi_hat <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    pi_hat[i] <- (1 + sum(p_values >= lambda)) / (n * (1 - lambda))
  }

  ambda_opt_index <- 1
  for (i in 2:length(pi_hat)) {
    if (pi_hat[i] > pi_hat[i - 1]) {
      lambda_opt_index <- i - 1
      break
    }
  }
  # Find the smallest lambda where pi_hat stops decreasing
  lambda_opt <- lambda_grid[lambda_opt_index]
  pi_hat_opt <- pi_hat[lambda_opt_index]

  # Step 2: Apply the BH procedure with the adjusted threshold
  p_sorted <- sort(p_values)
  k <- max(which(p_sorted <= (1:n) * q / (n * pi_hat_opt)))
  threshold <- p_sorted[k]

  # Return the results
  rejection_set <- which(p_values <= threshold)
  list(rejection_set = rejection_set, lambda_opt = lambda_opt, pi_hat_opt = pi_hat_opt)
}
