

HDMT_asy = function(p.alpha,p.beta,estws){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  nullprop <- estws
  c = nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  pi.01.est = nullprop$alpha01/c
  pi.10.est = nullprop$alpha10/c
  pi.00.est = nullprop$alpha00/c

  p.HDMT = pi.10.est*pmax+pi.01.est*pmax+pi.00.est*pmax^2

  p.HDMT[which(p.HDMT<=0)] = min(p.HDMT[which(p.HDMT>0)])
  p.HDMT[which(p.HDMT>=1)] = 1-1e-8

  return(p.HDMT)
}

