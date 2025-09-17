

HDMT_emp = function(p.alpha,p.beta,estws){
  input_pvalues = cbind(p.alpha,p.beta)
  pmax <- apply(input_pvalues,1,max)
  c = estws$alpha10+estws$alpha01+estws$alpha00
  pi.01.est = estws$alpha01/c
  pi.10.est = estws$alpha10/c
  pi.00.est = estws$alpha00/c

  x <- tryCatch(p_value_underH1(input_pvalues = input_pvalues, alpha1 = estws$alpha1, alpha2 = estws$alpha2), error = function(e) e)
  if(!inherits(x, "error")){
    correction = p_value_underH1(input_pvalues = input_pvalues, alpha1 = estws$alpha1, alpha2 = estws$alpha2)
    p_1j_H10 = correction[,1]
    p_2j_H01 = correction[,2]
    p.HDMT = pi.10.est*pmax*p_1j_H10+pi.01.est*pmax*p_2j_H01+pi.00.est*pmax^2
  }else if(inherits(x, "error")){
    p.HDMT = pi.10.est*pmax+pi.01.est*pmax+pi.00.est*pmax^2
  }

  p.HDMT[which(p.HDMT<=0)] = min(p.HDMT[which(p.HDMT>0)])
  p.HDMT[which(p.HDMT>=1)] = 1-1e-8

  return(p.HDMT)
}
