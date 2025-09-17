


##### functions
library(MASS)
#library(car)
library(locfdr)
library(fdrtool)
library(HDMT)
library(DACT)
library(lattice)
library(pracma)
`%!in%` <- Negate(`%in%`)

qqunif.plot<-function(pvalues,
                      should.thin=T, thin.obs.places=2, thin.exp.places=2,
                      xlab=expression(paste("Expected (",-log[10],~p[beta[j]],")")),
                      ylab=expression(paste("Observed (",-log[10],~p[beta[j]],")")),
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {


  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" ||
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }


  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }


  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }

  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()

  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }

  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points,
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}


## Sobel
Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
}

## MT-COMP
CompTestER = function(a,b,cdfX = NULL){
  ab = abs(a*b)
  var_vec = c(var(a),var(b))
  df_min = min(ab)/sqrt(max(var_vec))
  if(is.null(cdfX)){cdfX = CDF_p2(df_min, max(ab))}
  xout = ab/sqrt(sum(var_vec)-1)
  p_value = approx(x = cdfX$x, y = cdfX$y, xout = xout, yleft = 1, yright = 0, method = 'linear')$y
  return(list(pp=p_value,zz=safe_z(p_value)*sign(a*b)))
}

CDF_p2 = function(df_min, df_max){
  len_x = 1e5
  max_df_prod = df_max * 2
  grid_x = poly_space(0, max_df_prod, len_x, order = 10)
  grid_y = 2 * besselK(grid_x, 0) / pi
  integrand = (grid_y[1:(len_x - 1)] + grid_y[2:len_x])/2 * diff(grid_x)
  int_grid_y = rev(cumsum(rev(integrand))); int_grid_y[1] = 1
  cdfX = list(x = grid_x[1:(len_x - 1)], y = int_grid_y)
  return(cdfX)
}

safe_z=function(pp){ return(ifelse(pp<8e-324,40,qnorm(pp/2,lower.tail = F))) }

poly_space = function(a, b, n, order = 1){
  k = (b - a)^(1 - order)
  linsp = pracma::linspace(a, b, n)
  return(k * (linsp - a) ^ order + a)
}

## DACT origional
JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}

EfronCorrect = function(pval){
  z = stats::qnorm(1-pval)
  res <- locfdr(z,nulltype = 1)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}

nonnullPropEst <- function(x,u,sigma){
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) {

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

nullParaEst<-function (x,gamma=0.1){
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation

  n = length(x)
  t = c(1:1000)/200

  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0

  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)

  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }

  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]

  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)

  return(musigma=list(mu=uhat,s=shat))
}

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


#### HDMT
p_value_underH1 = function(input_pvalues, alpha1, alpha2, verbose = FALSE) {

  # input_pvalues contains two columns.
  # The first column is p_1j (p for alpha=0). The second column is p_2j (p for beta=0).
  # alpha1: proportion of null alpha=0
  # alpha2: proportion of null beta=0

  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)
  efdr1 <- rep(0,nmed)

  nmed  <- nrow(input_pvalues)
  cdf12 <- input_pvalues

  xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
  yy1 <- c(0,seq(1,nmed,by=1)/nmed)
  unique_xx1 <- unique(xx1)
  unique_yy1 <- yy1[which(!duplicated(xx1))]
  gfit1<- gcmlcm(unique_xx1,unique_yy1,type="lcm")
  xknots1 <- gfit1$x.knots[-1]
  Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)

  xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
  yy2 <- c(0,seq(1,nmed,by=1)/nmed)
  unique_xx2 <- unique(xx2)
  unique_yy2 <- yy2[which(!duplicated(xx2))]
  gfit2<- gcmlcm(unique_xx2,unique_yy2,type="lcm")
  xknots2 <- gfit2$x.knots[-1]
  Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)

  if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
  if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))


  orderq1 <- pmax
  orderq2 <- pmax

  gcdf1 <- pmax
  gcdf2 <- pmax
  for (i in 1:length(xknots1)) {
    if (i==1) {
      gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]]
    } else {
      if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
        if(verbose) {print(i)}
        temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]]
        gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
      }
    }
  }

  for (i in 1:length(xknots2)) {
    if (i==1) {
      gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]]
    } else {
      if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
        temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]]
        gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
      }
    }
  }


  gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
  gcdf2 <- ifelse(gcdf2>1,1,gcdf2)

  cdf12[,1] <- pmax(gcdf1,2) # p_1j under H10
  cdf12[,2] <- pmax(gcdf2,2) # p_2j under H01

  return(cdf12 = cdf12)

}

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

###correction
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

# Bisection method to find the largest t such that F_emp(t) <= significance_upper
# Define the search function using sorted index
search_sorted_index <- function(F_emp, Ts, control.method,significance_upper) {

  sorted_index <- order(Ts) # Get the indices that sort 'Ts'
  lower_bound <- 1
  upper_bound <- length(Ts)

  if(control.method == "size"){
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

    if(F_emp(t_s) - significance_upper < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }

  if(control.method == "FDR"){
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]

      if (F_emp(mid_value)*length(Ts)/mid_index < significance_upper ) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }

    t_s <- Ts[sorted_index[lower_bound]]

    if(F_emp(t_s)*length(Ts)/mid_index - significance_upper < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }

  }
}

alphas_estimate_tail <- function(X,estws,m0){
  X  <- as.matrix(X)
  input_pvalues = X
  nullprop <- estws

  m <- nrow(X)
  m0 <- ceiling(m0)
  r1 <- sort(X[, 1])[m0]
  r2 <- sort(X[, 2])[m0]

  solve.alpha <- function(x) {
    sum(log(X[X[, 1] <= r1, 1])) / m - min(0.99, nullprop$alpha1) * (r1*log(r1)-r1) - ((m0 / m  - r1 * nullprop$alpha1) / r1^{x} / (1 - nullprop$alpha1))*(1 - min(0.99, nullprop$alpha1)) * (r1^{x} * log(r1) - x^{-1} * r1^{x})
  }
  alpha1_hat <- uniroot(solve.alpha, c(1e-3, 1))$root

  solve.alpha <- function(x) {
    sum(log(X[X[, 2] <= r2, 2])) / m - min(0.99, nullprop$alpha2) * (r2*log(r2)-r2) - ((m0 / m  - r2 * nullprop$alpha2) / r2^{x} / (1 - nullprop$alpha2))*(1 - min(0.99, nullprop$alpha2)) * (r2^{x} * log(r2) - x^{-1} * r2^{x})
  }

  alpha2_hat <- uniroot(solve.alpha, c(1e-3, 1))$root

  C1 <- min((m0 / m  - r1 * nullprop$alpha1) / r1^{alpha1_hat} / (1 - nullprop$alpha1),3)
  C2 <- min((m0 / m  - r2 * nullprop$alpha2) / r2^{alpha2_hat} / (1 - nullprop$alpha2),3)

  return(c(alpha1_hat,alpha2_hat, C1, C2))
}


###correction
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


###correction

MDACT_pvalues = function(p.M,p.Y,estws){
  sim.num <- length(p.M)

  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-max(c,1-1e-6)

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

  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
  }

  pv <- vapply(Ts, Fnull, numeric(1))

  return(pv)
}

# Define the AS method based on the provided description
adaptive_bh_as <- function(p_values, q, delta = 0.01) {
  n <- length(p_values)

  # Step 1: Find the optimal lambda using the stopping rule
  lambda_grid <- seq(q, 1, by = delta)
  pi_hat <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    pi_hat[i] <- (1 + sum(p_values >= lambda)) / (n * (1 - lambda))
  }

  ambda_opt_index <- 1
  for (i in 2:length(pi_hat)) {
    if (pi_hat[i] > pi_hat[i - 1]) {
      lambda_opt_index <- i - 1
      break
    }
  }
  # Find the smallest lambda where pi_hat stops decreasing
  lambda_opt <- lambda_grid[lambda_opt_index]
  pi_hat_opt <- pi_hat[lambda_opt_index]

  # Step 2: Apply the BH procedure with the adjusted threshold
  p_sorted <- sort(p_values)
  k <- max(which(p_sorted <= (1:n) * q / (n * pi_hat_opt)))
  threshold <- p_sorted[k]

  # Return the results
  rejection_set <- which(p_values <= threshold)
  list(rejection_set = rejection_set, lambda_opt = lambda_opt, pi_hat_opt = pi_hat_opt)
}


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


index_fpn_fdr_power_func <- function(sat_index,true_set){
  if(length(true_set)==0){
    fpn <- ifelse(length(sat_index)>0,1,0)
    return(fpn)
  }else{
    fdr <-length(which(sat_index %!in% true_set))/max(length(sat_index),1)
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(fdr,power))
  }
}

size_func <- function(sat_index,true_set,sim.num){
  if(length(true_set)==0){
    fpn <- length(sat_index)
    return(fpn)
  }else{
    size <- length(which(sat_index %!in% true_set))/sim.num
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(size,power))
  }
}

generate_pvalues_indp_from_norm <- function(sim.num,ws,m,delta){
  pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)
  combinations <- as.matrix(expand.grid(rep(list(c(0, 1)),2)))[c(2,3,1,4),]
  # Calculate the number of units to be assigned to each set based on probabilities
  sim.num_per_set <- round(ws * sim.num)
  sets <- vector("list", length = length(ws))
  current_index <- 1
  for (i in 1:length(ws)) {
    if (sim.num_per_set[i] > 0) {
      sets[[i]] <- current_index:(current_index + sim.num_per_set[i] - 1)
      current_index <- current_index + sim.num_per_set[i]
    } else {
      sets[[i]] <- integer(0)  # Empty set
    }
  }
  mu_M <- sample(c(1,-1),length(c(unlist(sets[1]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),pras[1],1)
  mu_Y <- sample(c(1,-1),length(c(unlist(sets[2]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),pras[2],1)
  Z.M <- rnorm(sim.num,0,1); Z.M[c(unlist(sets[1]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),mu_M,1);
  Z.Y <- rnorm(sim.num,0,1); Z.Y[c(unlist(sets[2]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),mu_Y,1);
  p.M <- 2* (1 - pnorm(abs(Z.M))); p.Y <- 2* (1 - pnorm(abs(Z.Y)))
  p.M[p.M==0] <- 1e-17 ; p.Y[p.Y==0] <- 1e-17
  return(list(pvalues = cbind(p.M,p.Y),zvalues=cbind(Z.M,Z.Y)))
}


library(locfdr)
library(HDMT)
EstNull.func<-function (x,gamma=0.1){
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation

  n = length(x)
  t = c(1:1000)/200

  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0

  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)

  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }

  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]

  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)

  return(musigma=list(mu=uhat,s=shat))
}

epsest.func <- function(x,u,sigma){
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) {

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

JC_prop_estmate <- function(zvalues){
  musigma1=EstNull.func(zvalues[,1])
  musigma2=EstNull.func(zvalues[,2])
  pi1_first <- epsest.func(zvalues[,1],musigma1$mu,musigma1$s)
  pi2_first <- epsest.func(zvalues[,2],musigma2$mu,musigma2$s)
  return(c(pi1_first*(1-pi2_first),pi2_first*(1-pi1_first),(1-pi1_first)*(1-pi2_first),pi1_first*pi2_first))
}

# Function to reset temporary directory for each replication
reset_temp_dir <- function(tempdir_path) {
  # If the directory exists, delete it
  if (dir.exists(tempdir_path)) {
    unlink(tempdir_path, recursive = TRUE)  # Delete directory and its contents
  }

  # Recreate the temporary directory
  dir.create(tempdir_path, recursive = TRUE)

  # Set the directory for temporary files in R and bigstatsr
  Sys.setenv(TMPDIR = tempdir_path)
  options(bigstatsr.temporary_directory = tempdir_path)
}

