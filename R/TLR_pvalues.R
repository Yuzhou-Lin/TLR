TLR_pvalues = function(p.M,p.Y,estws,coefs,coefs2){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-min(c,0.5)

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M^(1-coefs[1]) * p.Y^(1-coefs[2]))

  ER00 <- function(t) {
    s <- max((t-pi.01.est*coefs2[1]^(-1)/coefs[1])/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)/coefs[1]+pi.10.est*coefs2[2]^(-1)),0)
    f <- function(x, pi.00.est, pi.01.est, pi.10.est, coefs,coefs2, t) {
      c <- coefs[1]^(1/(1-coefs[1]))*coefs[2]^(1/(1-coefs[2]))*(1/(1-coefs[2]))
      integrand <- ((t-pi.10.est*coefs2[2]^(-1)*x)/(pi.00.est*coefs2[1]^(-1)*coefs2[2]^(-1)*x+pi.01.est*coefs2[1]^(-1)))^(1/(1-coefs[1])) * x ^(coefs[2]/(1-coefs[2]))
      return(integrand*c)
    }
    up <- min((t/(pi.10.est*coefs2[2]^(-1))),1/coefs[2])
    re2 <- (s*coefs[2])^(1/(1-coefs[2]))
    if(s<=up){
      re1 <- integrate(f, lower = s, upper = up, pi.00.est = pi.00.est, pi.01.est = pi.01.est, pi.10.est = pi.10.est, coefs=coefs, coefs2=coefs2, t = t)
      return(re1$value + re2)
    }else{
      return(1)
    }
  }

  ER01 <- function(t) {
    s <- pmin(pmax((t - pi.10.est * coefs2[2]^(-1) * coefs[2]^{-1} * p.Y^{1 - coefs[2]}) /
                     (pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.Y^{1 - coefs[2]}), 0), 1)
    mean(s^{1 / (1 - coefs[1])})
  }

  ER10 <- function(t) {
    s <- pmin(pmax((t - pi.01.est * coefs2[1]^(-1) * coefs[1]^{-1} * p.M^{1 - coefs[1]}) /
                     (pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) +
                        pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * p.M^{1 - coefs[1]}), 0), 1)
    mean(s^{1 / (1 - coefs[2])})
  }

  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t)))
  }

  pv <- vapply(Ts, Fnull, numeric(1))

  return(pv)
}


