alphas_estimate_tail <- function(X,estws,m0){
  X  <- as.matrix(X)
  input_pvalues = X
  nullprop <- estws

  m <- nrow(X)
  m0 <- ceiling(m0)
  r1 <- sort(X[, 1])[m0]
  r2 <- sort(X[, 2])[m0]

  solve.alpha <- function(x) {
    sum(log(X[X[, 1] <= r1, 1])) / m - min(0.99, nullprop$alpha1) * (r1*log(r1)-r1) - ((m0 / m  - r1 * nullprop$alpha1) / r1^{x} / (1 - nullprop$alpha1))*(1 - min(0.99, nullprop$alpha1)) * (r1^{x} * log(r1) - x^{-1} * r1^{x})
  }
  alpha1_hat <- uniroot(solve.alpha, c(1e-3, 1))$root

  solve.alpha <- function(x) {
    sum(log(X[X[, 2] <= r2, 2])) / m - min(0.99, nullprop$alpha2) * (r2*log(r2)-r2) - ((m0 / m  - r2 * nullprop$alpha2) / r2^{x} / (1 - nullprop$alpha2))*(1 - min(0.99, nullprop$alpha2)) * (r2^{x} * log(r2) - x^{-1} * r2^{x})
  }

  alpha2_hat <- uniroot(solve.alpha, c(1e-3, 1))$root

  C1 <- min((m0 / m  - r1 * nullprop$alpha1) / r1^{alpha1_hat} / (1 - nullprop$alpha1),3)
  C2 <- min((m0 / m  - r2 * nullprop$alpha2) / r2^{alpha2_hat} / (1 - nullprop$alpha2),3)

  return(c(alpha1_hat,alpha2_hat, C1, C2))
}
