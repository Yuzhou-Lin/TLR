TLR_thre <- function(p.M, p.Y,coefs,coefs2,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)

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

  if(control.method == "FDR"){
    if(ifcond == 0){
      thr <- function(t) {
        c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                                pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                                (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                   pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    if(ifcond == 1){
      if(length(which(p.adjust(p.M,method = "fdr") < 0.1)) < 50 | length(which(p.adjust(p.Y,method = "fdr") < 0.1)) < 50){
        est10 <- 1
        est01 <- 1
        thr <- function(t) {
          (ER10(t) + ER01(t) - ER00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
        }
      }else{
        est10 <- mean(p.adjust(p.M,method = "fdr") < 0.1 & p.Y > 0.9)/punif(1-0.9)/mean(p.adjust(p.M,method = "fdr") < 0.1)
        est01 <- mean(p.adjust(p.Y,method = "fdr") < 0.1 & p.M > 0.9)/punif(1-0.9)/mean(p.adjust(p.Y,method = "fdr") < 0.1)
      }
      thr <- function(t) {
        (est10*ER10(t) + est01*ER01(t) +
           (pi.00.est - est10*(pi.00.est + pi.01.est) - est01*(pi.00.est + pi.10.est))*ER00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)

    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      if (thr(mid_value) < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]

    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }

  }
  if(control.method == "size"){
    thr <- function(t) {
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                          pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)

    F_emp <- function(t) {
      F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*ER10(t) +
                                   pi.01.est/(pi.01.est+pi.11.est)*ER01(t) +
                                   (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                      pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*ER00(t)))
      return(F_emp)
    }

    if(!inherits(x, "error")){
      t_s <- x
    }else{

      sorted_index <- order(Ts) # Get the indices that sort 'Ts'
      lower_bound <- 1
      upper_bound <- length(Ts)

      while (upper_bound - lower_bound > 1) {
        mid_index <- floor((lower_bound + upper_bound) / 2)
        mid_value <- Ts[sorted_index[mid_index]]

        if (F_emp(mid_value) < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }

    if(F_emp(round(t_s,6)) < significance_upper){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
}


