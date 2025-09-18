alphas_estimate_tail <- function(X,estws,m0){
  X  <- as.matrix(X)
  input_pvalues = X
  nullprop <- estws

  m <- nrow(X)
  m0 <- ceiling(m0)
  r1 <- sort(X[, 1])[m0]
  r2 <- sort(X[, 2])[m0]

  solve.alpha <- function(x) {
    sum(log(X[X[, 1] <= r1, 1])) / m - min(0.996, nullprop$alpha1) * (r1*log(r1)-r1) - ((m0 / m  - r1 * nullprop$alpha1) / r1^{x} / (1 - nullprop$alpha1))*(1 - min(0.996, nullprop$alpha1)) * (r1^{x} * log(r1) - x^{-1} * r1^{x})
  }

  alpha1_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if(inherits(alpha1_hat, "try-error")){
    x_vals <- seq(0.001, 1, by = 0.01)
    y_vals <- sapply(x_vals, solve.alpha)
    closest_index <- which.min(abs(y_vals))
    alpha1_hat <- x_vals[closest_index]
  }

  solve.alpha <- function(x) {
    sum(log(X[X[, 2] <= r2, 2])) / m - min(0.99, nullprop$alpha2) * (r2*log(r2)-r2) - ((m0 / m  - r2 * nullprop$alpha2) / r2^{x} / (1 - nullprop$alpha2))*(1 - min(0.99, nullprop$alpha2)) * (r2^{x} * log(r2) - x^{-1} * r2^{x})
  }

  alpha2_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if(inherits(alpha2_hat, "try-error")){
    x_vals <- seq(0.001, 1, by = 0.01)
    y_vals <- sapply(x_vals, solve.alpha)
    closest_index <- which.min(abs(y_vals))
    alpha2_hat <- x_vals[closest_index]
  }

  C1 <- max(0.1,min((m0 / m  - r1 * nullprop$alpha1) / r1^{alpha1_hat} / (1 - nullprop$alpha1),3))
  C2 <- max(0.1,min((m0 / m  - r2 * nullprop$alpha2) / r2^{alpha2_hat} / (1 - nullprop$alpha2),3))

  return(c(alpha1_hat,alpha2_hat, C1, C2))
}
