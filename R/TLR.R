TLR = function(p.M,p.Y,estws,coefs,coefs2,significance_upper,control.method,ifcond){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est = 1-c
  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c

  Ts <- pi.10.est * coefs2[2]^(-1) * coefs[2]^(-1) * (p.Y^(1-coefs[2])) +
    pi.01.est * coefs2[1]^(-1) * coefs[1]^(-1) * (p.M^(1-coefs[1])) +
    pi.00.est * coefs2[1]^(-1) * coefs2[2]^(-1) * coefs[1]^(-1) * coefs[2]^(-1) * (p.M^(1-coefs[1]) * p.Y^(1-coefs[2]))

  ss <- TLR_thre(p.M, p.Y,coefs,coefs2,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method,ifcond)

  return(ss)
}
