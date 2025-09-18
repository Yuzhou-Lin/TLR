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
