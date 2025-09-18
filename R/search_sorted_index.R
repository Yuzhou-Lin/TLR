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
