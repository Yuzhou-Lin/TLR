

poly_space = function(a, b, n, order = 1){
  k = (b - a)^(1 - order)
  linsp = pracma::linspace(a, b, n)
  return(k * (linsp - a) ^ order + a)
}
