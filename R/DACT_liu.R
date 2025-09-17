

DACT_liu <- function (p_a, p_b,Z_a= "NULL", Z_b= "NULL", correction = "NULL") {

  if(Z_a== "NULL" & Z_b == "NULL"){
    Z_a = stats::qnorm(p_a, lower.tail = F)
    Z_b = stats::qnorm(p_b, lower.tail = F)
  }

  pi0a = 1 - nonnullPropEst(Z_a, 0, 1)
  pi0b = 1 - nonnullPropEst(Z_b, 0, 1)
  if (pi0a > 1) {
    pi0a = 1
  }
  if (pi0b > 1) {
    pi0b = 1
  }
  p.mat = cbind(p_a, p_b)
  p3 = (apply(p.mat, 1, max))^2
  wg1 = pi0a * (1 - pi0b)
  wg2 = (1 - pi0a) * pi0b
  wg3 = pi0a * pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1, wg2, wg3)/wg.sum
  p_dact = wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
  if (correction == "Efron") {
    x <- tryCatch(EfronCorrect(p_dact), error = function(e) e)
    if(inherits(x, "error")){
      p_dact = NULL
    }else{
      p_dact = EfronCorrect(p_dact)
    }
  }
  if (correction == "JC") {
    x <- tryCatch(JCCorrect(p_dact), error = function(e) e)
    if(inherits(x, "error")){
      p_dact = NULL
    }else{
      p_dact = JCCorrect(p_dact)
    }
  }
  return(p_dact)
}
