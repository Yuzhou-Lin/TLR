MDACT = function(p.M,p.Y,estws,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est

  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c

  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2

  ss <- DACT_thr(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond)

  return(ss)
}


