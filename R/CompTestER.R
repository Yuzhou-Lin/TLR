

CompTestER = function(a,b,cdfX = NULL){
  ab = abs(a*b)
  var_vec = c(var(a),var(b))
  df_min = min(ab)/sqrt(max(var_vec))
  if(is.null(cdfX)){cdfX = CDF_p2(df_min, max(ab))}
  xout = ab/sqrt(sum(var_vec)-1)
  p_value = approx(x = cdfX$x, y = cdfX$y, xout = xout, yleft = 1, yright = 0, method = 'linear')$y
  return(list(pp=p_value,zz=safe_z(p_value)*sign(a*b)))
}

