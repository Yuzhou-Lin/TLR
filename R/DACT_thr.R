DACT_thr <- function(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)

  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2

  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }

    x <- tryCatch(integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- integrate(F_0int, 0, 1)
    }
    res$value
  }

  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }

  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }

  if(control.method == "FDR"){
    if(ifcond == 0){
      thr <- function(t) {
        c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                                pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                                (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                   pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    if(ifcond == 1){
      if(length(which(p.adjust(p.M,method = "fdr") < 0.1)) < 50 | length(which(p.adjust(p.Y,method = "fdr") < 0.1)) < 50){
        est10 <- 1
        est01 <- 1
        thr <- function(t) {
          (F10(t) + F01(t) - F00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
        }
      }else{
        est10 <- mean(p.adjust(p.M,method = "fdr") < 0.1 & p.Y > 0.9)/punif(1-0.9)/mean(p.adjust(p.M,method = "fdr") < 0.1)
        est01 <- mean(p.adjust(p.Y,method = "fdr") < 0.1 & p.M > 0.9)/punif(1-0.9)/mean(p.adjust(p.Y,method = "fdr") < 0.1)
      }
      thr <- function(t) {
        (est10*F10(t) + est01*F01(t) +
           (pi.00.est - est10*(pi.00.est + pi.01.est) - est01*(pi.00.est + pi.10.est))*F00(t)) / max(mean(Ts < t), 1 / sim.num) - significance_upper
      }
    }
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    thr1 <- 50
    while (upper_bound - lower_bound > thr1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error")){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound - (thr1/10*3)]]

    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }

  }
  if(control.method == "size"){
    thr <- function(t) {
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                          pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)
    F_emp <- function(t) {
      F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                                   pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                                   (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                                      pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
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

        x <- tryCatch(F_emp(round(mid_value,6)), error = function(e) e)

        if(!inherits(x, "error")){
          fn <- F_emp(round(mid_value,6))
        }else{
          fn <- 0
        }

        if (fn < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }

    x <- tryCatch(F_emp(t_s), error = function(e) e)
    if(!inherits(x, "error")){
      if(F_emp(t_s) < significance_upper){
        return(which(Ts <= t_s))
      }else{
        return(which(Ts < t_s))
      }
    }else{
      return(which(Ts <= t_s))
    }
  }
}
