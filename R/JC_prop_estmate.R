

JC_prop_estmate <- function(zvalues){
  musigma1=EstNull.func(zvalues[,1])
  musigma2=EstNull.func(zvalues[,2])
  pi1_first <- epsest.func(zvalues[,1],musigma1$mu,musigma1$s)
  pi2_first <- epsest.func(zvalues[,2],musigma2$mu,musigma2$s)
  return(c(pi1_first*(1-pi2_first),pi2_first*(1-pi1_first),(1-pi1_first)*(1-pi2_first),pi1_first*pi2_first))
}
