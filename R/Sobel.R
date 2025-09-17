


## Sobel
Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
}
