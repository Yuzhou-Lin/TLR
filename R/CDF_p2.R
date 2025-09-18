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

