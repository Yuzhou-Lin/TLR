

HDMT_control = function(p.M,p.Y,estws,coefs, exact = 1,significance_upper, control.method){
  sim.num <- length(p.M)

  F_emp <- function(x){x}

  if(exact==1){
    Ts <- HDMT_emp(p.M,p.Y,estws)
  }
  if(exact==0){
    Ts <- HDMT_asy(p.M,p.Y,estws)
  }

  ss <- search_sorted_index(F_emp,Ts,control.method,significance_upper)

  return(ss)
}
